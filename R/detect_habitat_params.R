#' Derive DBSCAN habitat-delineation parameters from the data
#'
#' Estimates the two DBSCAN parameters (`eps` and `minPts`) used by [dbscan_habitats()]
#' from empirical tracking data, removing the need for hard-coded assumptions. Both parameters
#' are derived at the population level and have an explicit behavioural meaning. This
#' generalized version dynamically adapts to whether behavioural state tracking data is present
#' or absent.
#'
#' @section eps -- the within-site movement scale:
#' `eps` represents the radius linking telemetry fixes that belong to the same temporary
#' habitat patch. It can be derived in two modes:
#' \itemize{
#'   \item **State-Filtered Mode (default):** When `track` carries the `state_col` column
#'   (e.g. BALM `cluster_type`) and `stationary_states` is set, consecutive great-circle
#'   steps whose *both* endpoints are stationary (e.g. BALM "LL") are isolated. This cleanly
#'   strips away high-speed flights and between-site transits before fitting.
#'   \item **Generalized Mode (automatic fallback):** If `state_col` is absent from `track`
#'   (or no `stationary_states` are given), all valid positive step lengths across the
#'   population are pooled instead. This lets non-migratory species be analysed directly on
#'   the raw track, without a prior BALM segmentation step.
#' }
#' The function defaults to State-Filtered Mode but falls back to Generalized Mode
#' automatically (with a message) when the state column is not present, so the same call
#' works for both BALM-segmented and raw tracks.
#'
#' In both modes, a Gaussian Mixture Model (GMM) is fitted to the log-transformed step lengths.
#' The mixture is partitioned at the **largest gap** between consecutive component means,
#' splitting the distribution into a lower "within-site complex" (resting, foraging, local transits)
#' and higher relocation modes. `eps` is set as the `eps_quantile` (default 0.95) percentile
#' of this within-site complex. Pooling across the population ensures statistical stability and
#' prevents the overfitting (spurious components) common in per-individual fits.
#'
#' Note: in Generalized Mode the within-site/relocation split relies entirely on the GMM
#' largest-gap heuristic. This is robust for sedentary species, but for central-place
#' foragers (residents that commute between a roost and feeding sites) the commute steps can
#' inflate `eps`; consider a light stationarity filter, or BALM, in that case.
#'
#' @section minPts -- the minimum-residence floor:
#' `minPts` defines the minimum number of local fixes a cluster must contain to be recognized
#' as a valid habitat. By default, it is determined via the Kneedle geometry algorithm, which
#' captures the "knee" of the aggregate habitat-count vs. `minPts` curve at the derived `eps`.
#' This curve marks the physical transition from culling accidental noise points to eroding
#' genuine habitat clusters.
#'
#' Its operational meaning is temporal: `minPts * sampling_interval = minimum residence duration`.
#' Alternatively, users can explicitly supply `min_residence` (in days) to enforce a direct
#' behavioral floor, which is converted to fix counts using the median sampling interval.
#'
#' @param x An `sf` POINT object containing the stationary subset of points to be clustered
#'   (the data over which `minPts` is scanned).
#' @param track An `sf` POINT object containing the complete telemetry trajectory. Used to
#'   extract step lengths, calculate sampling intervals, and analyze movement states. Required
#'   unless a custom `eps` value is explicitly provided without specifying `min_residence`.
#' @param individual_col Character string identifying the column for unique individual IDs
#'   in `x` and `track`. If `NULL`, all rows are treated as a single group.
#' @param time_col Character string specifying the timestamp column in `track`.
#' @param state_col Movement-state label column in `track` (default `"cluster_type"`, the BALM
#'   output). If this column is **absent** from `track`, state filtering is bypassed and the
#'   function falls back to Generalized Mode (all steps pooled) -- so raw, un-segmented tracks
#'   from non-migratory species work without change. Set `state_col = NULL` to force the
#'   generalized mode even when a state column exists.
#' @param stationary_states Vector of state labels treated as local/stationary for the `eps`
#'   step set (default `c("LL")`). If `NULL`, state-based filtering is bypassed.
#' @param flight_states Vector of state labels treated as flight/relocation (default `c("HH")`),
#'   used only for the regime diagnostic (`gap_ratio`).
#' @param eps Optional numeric value. If supplied, the empirical derivation of `eps` (and
#'   the tracking data requirements tied to it) is skipped.
#' @param eps_quantile Numeric percentile (0 to 1) extracted from the within-site complex
#'   to define `eps` (default `0.95`).
#' @param minpts_max Integer establishing the upper limit for the `minPts` knee scan (default `60`).
#' @param min_residence Optional minimum patch residency duration expressed in **days**. If
#'   supplied, `minPts` bypasses the knee calculation and directly computes fix thresholds.
#' @param verbose Logical; if `TRUE`, outputs the derived values and diagnostic logs to the
#'   console (default `TRUE`).
#'
#' @return A named list containing:
#' \itemize{
#'   \item `eps`: Calculated or user-supplied neighborhood distance threshold (km).
#'   \item `minPts`: Calculated or converted minimum point density threshold (integer).
#'   \item `regime`: Diagnostic label indicating behavior pooling (`"migratory"`,
#'         `"continuous"`, `"all-data-pooled"`, or `"user-supplied"`).
#'   \item `gap_ratio`: The ratio of median flight step lengths to `eps` (diagnostic asset).
#'   \item `n_components`: Number of GMM components selected by BIC during clustering.
#'   \item `sampling_interval_h`: Median telemetry sampling interval expressed in hours.
#'   \item `min_residence_days`: The minimum habitat residence duration implied by `minPts`.
#' }
#'
#' @seealso [dbscan_habitats()].
#' @examples
#' \dontrun{
#' # Example 1: BALM-segmented track -- State-Filtered Mode is used by default
#' p1 <- detect_habitat_params(stat_sf, track = seg,
#'                             individual_col = "Individual", time_col = "Time")
#'
#' # Example 2: Raw track from a non-migratory species (no state column) --
#' # the function falls back to Generalized Mode automatically
#' p2 <- detect_habitat_params(raw_sf, track = raw_sf,
#'                             individual_col = "Individual", time_col = "Time")
#' }
#' @importFrom mclust Mclust mclustBIC
#' @export
detect_habitat_params <- function(x, track = NULL, individual_col = NULL,
                                  time_col = NULL, state_col = "cluster_type",
                                  stationary_states = c("LL"),
                                  flight_states = c("HH"),
                                  eps = NULL, eps_quantile = 0.95,
                                  minpts_max = 60L, min_residence = NULL,
                                  verbose = TRUE) {
  
  # ---- Input Validation ----
  if (!inherits(x, "sf"))
    stop("`x` must be an sf object.", call. = FALSE)
  
  for (pk in c("dbscan", "sf")) {
    if (!requireNamespace(pk, quietly = TRUE))
      stop(sprintf("Package '%s' is required.", pk), call. = FALSE)
  }
  
  rad <- pi / 180
  regime <- "user-supplied"
  gap_ratio <- NA_real_
  n_comp <- NA_integer_
  dt_med <- NA_real_
  
  # ---- Track Demands & Analysis Flow ----
  need_track <- is.null(eps) || !is.null(min_residence)
  if (need_track) {
    if (is.null(track) || !inherits(track, "sf"))
      stop("`track` (the trajectory) is required to derive `eps` or to honour `min_residence`.", call. = FALSE)
    if (is.null(time_col) || !time_col %in% names(track))
      stop("`time_col` must name a valid column in `track`.", call. = FALSE)
    
    tgrp <- if (is.null(individual_col)) rep(1L, nrow(track)) else track[[individual_col]]
    tco  <- sf::st_coordinates(track)
    tt   <- track[[time_col]]
    
    # Decide whether behavioural state-filtering applies. The default `state_col`
    # is the BALM column; if it is absent from `track` we silently fall back to
    # Generalized Mode, so raw (un-segmented) tracks work with the same call.
    use_states <- !is.null(state_col) && (state_col %in% names(track)) && !is.null(stationary_states)
    if (verbose && !use_states) {
      if (!is.null(state_col) && !(state_col %in% names(track)))
        message(sprintf(
          "[detect_habitat_params] state column '%s' not found in `track`; falling back to generalized (all-step) mode.",
          state_col))
      else
        message("[detect_habitat_params] no state filtering; using generalized (all-step) mode.")
    }
    if (use_states) {
      tst <- as.character(track[[state_col]])
    }
    
    loc_steps <- numeric(0)
    fl_steps  <- numeric(0)
    dts       <- numeric(0)
    
    # Process steps per individual
    for (gg in unique(tgrp)) {
      sel <- which(tgrp == gg)
      if (length(sel) < 2L) next
      
      # Enforce chronological sorting
      o   <- sel[order(tt[sel])]
      s   <- .gc_steps(tco[o, , drop = FALSE])
      dts <- c(dts, as.numeric(diff(tt[o]), units = "hours"))
      if (!length(s)) next
      
      if (use_states) {
        # State-Filtered Mode: keep steps whose BOTH endpoints are stationary
        st <- tst[o]
        s1 <- st[-length(st)]
        s2 <- st[-1]
        kl <- s1 %in% stationary_states & s2 %in% stationary_states
        loc_steps <- c(loc_steps, s[kl & s > 0])
        
        if (!is.null(flight_states)) {
          kf <- s1 %in% flight_states | s2 %in% flight_states
          fl_steps  <- c(fl_steps,  s[kf & s > 0])
        }
      } else {
        # Generalized Mode: accumulate all valid positive steps
        loc_steps <- c(loc_steps, s[s > 0])
      }
    }
    
    # Extract the population-level median sampling interval
    dt_med <- stats::median(dts[dts > 0], na.rm = TRUE)
    
    # ---- Mathematical Derivation of eps ----
    if (is.null(eps)) {
      if (!requireNamespace("mclust", quietly = TRUE))
        stop("Package 'mclust' is required to derive `eps`.", call. = FALSE)
      
      er  <- .habitat_eps_local(loc_steps, eps_quantile)
      eps <- er$eps
      n_comp <- er$n_components
      
      if (is.na(eps))
        stop("Could not derive `eps` (too few valid tracking steps available).", call. = FALSE)
      
      # Assign movement classification tags based on setup
      if (use_states && length(fl_steps) > 0L) {
        regime    <- "migratory"
        gap_ratio <- stats::median(fl_steps) / eps
      } else {
        regime    <- if (use_states) "continuous" else "all-data-pooled"
        gap_ratio <- NA_real_
      }
    }
  }
  
  # ---- Derivation of minPts ----
  grp    <- if (is.null(individual_col)) rep(1L, nrow(x)) else x[[individual_col]]
  ids    <- unique(grp)
  co_all <- sf::st_coordinates(x)
  
  if (!is.null(min_residence)) {
    # Method A: Direct temporal conversion based on user inputs
    if (is.na(dt_med) || dt_med <= 0)
      stop("Cannot convert `min_residence` without a valid tracking sampling interval.", call. = FALSE)
    minPts <- max(2L, as.integer(round(min_residence * 24 / dt_med)))
  } else {
    # Method B: Kneedle geometric optimization across standard thresholds
    mg  <- 2:as.integer(minpts_max)
    agg <- vapply(mg, function(m) {
      sum(vapply(ids, function(gg) {
        sel <- which(grp == gg)
        if (length(sel) < m) return(0L)
        co   <- co_all[sel, , drop = FALSE]
        latm <- mean(co[, 2])
        # Project local coordinates locally to approximate spatial distance metrics
        xy   <- cbind(co[, 1] * cos(latm * rad) * 111.320, co[, 2] * 111.320)
        cl   <- dbscan::dbscan(xy, eps = eps, minPts = m)$cluster
        length(unique(cl[cl > 0]))
      }, integer(1)))
    }, integer(1))
    
    minPts <- .knee(mg, agg)
    if (is.na(minPts))
      stop("Could not derive `minPts` (degenerate/flat habitat-count curve).", call. = FALSE)
  }
  
  res_days <- if (!is.na(dt_med) && dt_med > 0) minPts * dt_med / 24 else NA_real_
  
  # ---- Console Reporting Output ----
  if (verbose) {
    if (identical(regime, "user-supplied")) {
      message(sprintf("eps    = %.3f km  (user-supplied)", eps))
    } else {
      gmsg <- if (!is.na(gap_ratio)) sprintf("flight/eps gap %.0fx", gap_ratio)
      else if (identical(regime, "all-data-pooled")) "generalized: all steps pooled"
      else "states used, no flight class"
      message(sprintf(
        "eps    = %.3f km  (%s; %.2f-quantile of %d-comp GMM; %s)",
        eps, regime, eps_quantile, n_comp, gmsg))
    }
    
    if (!is.na(res_days)) {
      message(sprintf("minPts = %d  (~%.2f d min residence at %.1f h sampling; %s)",
                      minPts, res_days, dt_med,
                      if (!is.null(min_residence)) "user min_residence" else "knee"))
    } else {
      message(sprintf("minPts = %d", minPts))
    }
  }
  
  # Return final operational configurations
  list(eps = eps, minPts = minPts, regime = regime, gap_ratio = gap_ratio,
       n_components = n_comp, sampling_interval_h = dt_med,
       min_residence_days = res_days)
}
# ---------------------------------------------------------------------------
# Internal Helper Functions
# ---------------------------------------------------------------------------
#' Great-circle step lengths (km) for time-ordered coordinates
#' @keywords internal
#' @noRd
.gc_steps <- function(co) {
  n <- nrow(co)
  if (n < 2L) return(numeric(0))
  lon <- co[, 1]
  lat <- co[, 2]
  rad <- pi / 180
  dlat <- diff(lat) * rad
  dlon <- diff(lon) * rad
  a <- sin(dlat / 2)^2 + cos(lat[-n] * rad) * cos(lat[-1] * rad) * sin(dlon / 2)^2
  6371 * 2 * asin(pmin(1, sqrt(a)))
}
#' Derive eps from pooled data using a localized Log-Gaussian Mixture Model split
#' @keywords internal
#' @noRd
.habitat_eps_local <- function(steps, q = 0.95) {
  steps <- steps[is.finite(steps) & steps > 0]
  if (length(steps) < 30L)
    return(list(eps = NA_real_, n_components = NA_integer_))
  
  # Sort to ensure invariant inputs against initial configuration shifts
  lv <- sort(log(steps))
  mclustBIC <- mclust::mclustBIC
  
  # Preserve original environment RNG states while keeping GMM setups deterministic
  .old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(.old_seed))
    assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(1L)
  
  fit <- tryCatch(mclust::Mclust(lv, G = 1:5, verbose = FALSE),
                  error = function(e) NULL)
  
  if (is.null(fit))
    return(list(eps = as.numeric(exp(stats::quantile(lv, q, names = FALSE))),
                n_components = 1L))
  
  G  <- fit$G
  mu <- fit$parameters$mean
  v  <- fit$parameters$variance$sigmasq
  if (length(v) == 1L) v <- rep(v, G)
  w  <- fit$parameters$pro
  
  # Order components by ascending mean
  om <- order(mu); mu <- mu[om]; v <- v[om]; w <- w[om]
  
  if (G == 1L)
    return(list(eps = as.numeric(exp(stats::quantile(lv, q, names = FALSE))),
                n_components = 1L))
  
  # Identify separation index targeting the largest component gap
  cut <- which.max(diff(mu))
  loc <- seq_len(cut)
  wl  <- w[loc] / sum(w[loc])
  
  xs  <- seq(min(lv), mu[cut + 1L], length.out = 4000)
  cdf <- vapply(xs, function(z) sum(wl * stats::pnorm(z, mu[loc], sqrt(v[loc]))),
                numeric(1))
  
  list(eps = exp(xs[which.min(abs(cdf - q))]), n_components = G)
}
#' Kneedle geometric estimation locating the furthest orthogonal value from the chord vector
#' @keywords internal
#' @noRd
.knee <- function(x, y) {
  if (length(x) < 3L || diff(range(y)) == 0) return(NA_integer_)
  
  # Scale values to a normalized unit square
  xs <- (x - min(x)) / (max(x) - min(x))
  ys <- (y - min(y)) / diff(range(y))
  
  x1 <- xs[1]; y1 <- ys[1]
  x2 <- xs[length(xs)]; y2 <- ys[length(ys)]
  
  # Perpendicular distance from each point to the end-to-end chord
  d  <- abs((y2 - y1) * xs - (x2 - x1) * ys + x2 * y1 - y2 * x1) /
    sqrt((y2 - y1)^2 + (x2 - x1)^2)
  
  as.integer(x[which.max(d)])
}