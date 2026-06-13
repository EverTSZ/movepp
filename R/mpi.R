#' Detect the Dominant Migration Axis via PCA
#'
#' Automatically determines whether a tracking dataset is dominated by
#' latitudinal, longitudinal, or mixed (diagonal) movement, eliminating
#' the need for users to know this property *a priori*. Principal
#' Components Analysis (PCA) is applied to the (lon, lat) coordinates;
#' the fraction of total variance captured by PC1 determines whether the
#' movement is treated as unidirectional or mixed.
#'
#' @param points An `sf` POINT object.
#' @param pc1_threshold Numeric in (0, 1]; fraction of variance PC1 must
#'   explain to call the movement unidirectional (default 0.85).
#'
#' @return A list of class `"movepp_dominant_axis"` with elements
#'   `primary_axis`, `secondary_axis` (or `NULL`),
#'   `pc1_variance_explained`, `pc1_loading_lon`, `pc1_loading_lat`,
#'   and the underlying `pca` (`prcomp`) object.
#'
#' @examples
#' \dontrun{
#' demo <- make_demo_track()
#' detect_dominant_axis(demo)
#' }
#'
#' @importFrom stats prcomp complete.cases
#' @export
detect_dominant_axis <- function(points, pc1_threshold = 0.85) {
  if (!inherits(points, "sf")) {
    stop("`points` must be an `sf` object.", call. = FALSE)
  }
  if (pc1_threshold <= 0 || pc1_threshold > 1) {
    stop("`pc1_threshold` must be in (0, 1].", call. = FALSE)
  }
  
  coords <- sf::st_coordinates(points)
  coords <- coords[stats::complete.cases(coords), , drop = FALSE]
  
  if (nrow(coords) < 3L) {
    stop("Need at least 3 points to compute PCA.", call. = FALSE)
  }
  
  pca <- stats::prcomp(coords, scale. = FALSE, center = TRUE)
  
  pc1_var <- pca$sdev[1]^2 / sum(pca$sdev^2)
  pc1_lon <- abs(pca$rotation["X", "PC1"])
  pc1_lat <- abs(pca$rotation["Y", "PC1"])
  
  result <- list(
    pc1_variance_explained = pc1_var,
    pc1_loading_lon        = pc1_lon,
    pc1_loading_lat        = pc1_lat,
    pca                    = pca
  )
  
  if (pc1_var >= pc1_threshold) {
    if (pc1_lat > pc1_lon) {
      result$primary_axis   <- "latitude"
      result$secondary_axis <- NULL
    } else {
      result$primary_axis   <- "longitude"
      result$secondary_axis <- NULL
    }
  } else {
    if (pc1_lat > pc1_lon) {
      result$primary_axis   <- "latitude"
      result$secondary_axis <- "longitude"
    } else {
      result$primary_axis   <- "longitude"
      result$secondary_axis <- "latitude"
    }
  }
  
  class(result) <- "movepp_dominant_axis"
  result
}


#' Print method for movepp_dominant_axis objects
#'
#' @param x A `movepp_dominant_axis` object from [detect_dominant_axis()].
#' @param ... Ignored.
#' @export
print.movepp_dominant_axis <- function(x, ...) {
  cat("Dominant Movement Axis (PCA-based)\n")
  cat("-----------------------------------\n")
  cat(sprintf("Primary axis:   %s\n", x$primary_axis))
  cat(sprintf("Secondary axis: %s\n",
              if (is.null(x$secondary_axis)) "(none, unidirectional)"
              else x$secondary_axis))
  cat(sprintf("PC1 explains:   %.1f%% of variance\n",
              100 * x$pc1_variance_explained))
  cat(sprintf("PC1 loadings:   lon = %.3f, lat = %.3f\n",
              x$pc1_loading_lon, x$pc1_loading_lat))
  invisible(x)
}


#' Compute the Migration Phase Index (MPI)
#'
#' Calculates a per-record Migration Phase Index that combines *where* a
#' fix sits along the population's dominant migration axis with *when* in
#' the seasonal cycle it occurred. MPI is a scalar in \[0, 1] that
#' approaches 1 during the breeding residency at the breeding end of the
#' axis and approaches 0 during the wintering residency at the opposite
#' end. It is fully data-driven: nothing about the species' geometry,
#' hemisphere, or calendar is hard-coded.
#'
#' @section Spatial term (per individual):
#' The migration axis is found by PCA on `(lon, lat)`. If PC1 explains at
#' least `pc1_threshold` of the variance the movement is treated as
#' unidirectional (primary axis only); otherwise PC2 is also used
#' (mixed-axis migration). Each individual is normalised against *its
#' own* axis extremes, so partial migrants and individuals with
#' different geographic ranges are placed on a common 0--1 scale. This is
#' the "space = individual" anchor.
#'
#' The extremes used for normalisation are the **near-maximal** quantiles
#' `anchor_q` and `1 - anchor_q` (default 0.02, i.e. the 2nd and 98th
#' percentiles), *not* a moderate tail like the 10th/90th. The breeding
#' extreme is a *spatial* limit that can be badly under-represented in
#' *time*: an individual that skips breeding in some years spends most of
#' its northern time at a lower, non-breeding settlement, which would drag
#' a 90th-percentile anchor down and clamp the true breeding extreme to 1
#' -- making a non-breeding stop indistinguishable from real breeding.
#' Near-maximal anchors reserve `p_norm = 1` for the genuine extreme. The
#' separate, more moderate `q` tail is used only to collect the fixes
#' "at the extreme" whose dates set the population temporal anchors below.
#'
#' @section Temporal anchors (population level):
#' Which end of the axis is the *breeding* end is decided by **temporal
#' concentration**, not by latitude or photoperiod (both have too many
#' counterexamples across migration geometries). For each axis extreme,
#' the day-of-year of every fix near that extreme is pooled across the
#' population and its circular concentration `R` (mean resultant length)
#' is computed. Breeding sites are occupied in a tight seasonal window,
#' wintering sites over a longer, more diffuse period, so the extreme
#' with the **higher `R`** is the breeding end. A poleward cross-check is
#' used only to break near-ties. The breeding and wintering day-of-year
#' anchors (`breed_doy`, `winter_doy`) are the circular means of the
#' pooled fix dates at each end; they may be overridden directly. These
#' are the "time = population" anchors, derived from the same fix-sets
#' that define the spatial extremes, so space and time stay consistent.
#'
#' @section Temporal term (triangular tent):
#' The seasonal term `d_norm` is a triangular function of day-of-year,
#' anchored at the two population day-of-year points: it is 0 at the
#' wintering anchor, rises linearly to 1 at the breeding anchor, then
#' falls linearly back to 0. The two arcs (spring up, autumn down) fold
#' onto the same value, so MPI is a direction-agnostic breeding-phase
#' index -- spring and autumn migration receive the same score by design.
#'
#' @section MPI:
#' Single-axis: \eqn{\mathrm{MPI} = \sqrt{p_{\mathrm{norm}}\,
#' d_{\mathrm{norm}}}}. Mixed-axis:
#' \eqn{\mathrm{MPI} = \sqrt[3]{p_{\mathrm{norm}}\, s_{\mathrm{norm}}\,
#' d_{\mathrm{norm}}}}.
#'
#' @param points An `sf` POINT object (typically habitat-labelled
#'   stationary fixes).
#' @param individual_col Column identifying individuals. Spatial
#'   normalisation is per individual; temporal anchors are pooled across
#'   individuals. If `NULL`, all rows are treated as one individual.
#' @param time_col Column giving fix times (POSIXct/POSIXt/Date).
#'   Required.
#' @param q Tail fraction defining which fixes count as "at the breeding
#'   / wintering extreme" for pooling occupancy dates (the population
#'   temporal anchors and polarity test). Default 0.10. This does **not**
#'   set the spatial normalisation anchors -- see `anchor_q`.
#' @param anchor_q Near-maximal tail fraction for the per-individual
#'   spatial normalisation anchors (default 0.02; `p_norm = 0` at the
#'   `anchor_q` quantile, `1` at the `1 - anchor_q` quantile). Smaller =
#'   anchors closer to the true extreme (captures intermittent breeding
#'   but more sensitive to outlier fixes).
#' @param min_span_frac An individual contributes to the population
#'   temporal anchors only if its span along the axis is at least this
#'   fraction of the population span (default 0.30). Resident / barely
#'   moving individuals are thus excluded from anchor voting but still
#'   receive an MPI.
#' @param winter_doy,breed_doy Optional day-of-year anchors (1--365). If
#'   supplied, the corresponding date is used directly for the temporal tent
#'   **and** to orient the spatial polarity (which axis extreme is the
#'   breeding end is set to the extreme whose occupancy dates sit nearest
#'   `breed_doy` / farthest from `winter_doy`), so the spatial and temporal
#'   terms stay consistent. This overrides the automatic
#'   temporal-concentration polarity test -- use it whenever the printed
#'   polarity disagrees with the species' phenology.
#' @param pc1_threshold Variance fraction above which movement is treated
#'   as unidirectional (default 0.85).
#' @param verbose Logical; print the derived axis, anchors and polarity.
#'
#' @return The input `sf` with an added numeric `mpi` column in \[0, 1],
#'   plus attributes `mpi_winter_doy`, `mpi_breed_doy` (day-of-year) and
#'   `mpi_pc1_var`.
#'
#' @examples
#' \dontrun{
#' hab_mpi <- compute_mpi(hab, individual_col = "Individual",
#'                        time_col = "Time")
#' hist(hab_mpi$mpi)
#' }
#'
#' @seealso [detect_dominant_axis()] for the underlying axis detection.
#'
#' @importFrom stats prcomp quantile complete.cases
#' @export
compute_mpi <- function(points,
                        individual_col = NULL,
                        time_col       = NULL,
                        q              = 0.10,
                        anchor_q       = 0.02,
                        min_span_frac  = 0.30,
                        winter_doy     = NULL,
                        breed_doy      = NULL,
                        pc1_threshold  = 0.85,
                        verbose        = TRUE) {
  
  if (!inherits(points, "sf"))
    stop("`points` must be an `sf` object.", call. = FALSE)
  if (is.null(time_col) || !time_col %in% names(points))
    stop("`time_col` is required and must name a column in `points`.",
         call. = FALSE)
  tvec <- points[[time_col]]
  if (!inherits(tvec, c("POSIXct", "POSIXt", "Date")))
    stop("`time_col` must reference a POSIXct/POSIXt/Date column.",
         call. = FALSE)
  
  doy <- as.integer(format(tvec, "%j")) %% 365L      # 0..364 (366 -> 1)
  co  <- sf::st_coordinates(points)
  lon <- co[, "X"]; lat <- co[, "Y"]
  grp <- if (is.null(individual_col)) rep(1L, nrow(points)) else
    points[[individual_col]]
  
  # ---- circular helpers (day-of-year, period 365) ----------------------
  rad <- 2 * pi / 365
  circ_mean <- function(d) {
    d <- d[is.finite(d)]
    if (!length(d)) return(NA_real_)
    a <- d * rad
    (atan2(mean(sin(a)), mean(cos(a))) / rad) %% 365
  }
  circ_R <- function(d) {
    d <- d[is.finite(d)]
    if (length(d) < 3L) return(NA_real_)
    a <- d * rad
    sqrt(mean(cos(a))^2 + mean(sin(a))^2)
  }
  
  # ---- PCA axes (population) -------------------------------------------
  ok  <- stats::complete.cases(lon, lat)
  if (sum(ok) < 3L)
    stop("Fewer than 3 complete coordinates; cannot derive a migration ",
         "axis.", call. = FALSE)
  pca <- stats::prcomp(cbind(lon, lat)[ok, , drop = FALSE],
                       center = TRUE, scale. = FALSE)
  pc1_var <- pca$sdev[1]^2 / sum(pca$sdev^2)
  ctr     <- pca$center
  scores  <- cbind(lon - ctr[1], lat - ctr[2]) %*% pca$rotation
  ax1 <- scores[, 1]
  use_secondary <- pc1_var < pc1_threshold
  ax2 <- if (use_secondary) scores[, 2] else NULL
  
  # ---- per-axis processing ---------------------------------------------
  # Returns per-row p in [0,1] (1 = breeding end) plus the pooled
  # breeding/wintering fix dates used to set the temporal anchors.
  process_axis <- function(ax) {
    p_norm  <- rep(NA_real_, length(ax))
    hi_dates <- numeric(0); lo_dates <- numeric(0)
    hi_lat   <- numeric(0); lo_lat   <- numeric(0)
    pop_span <- diff(range(ax[is.finite(ax)]))
    if (!is.finite(pop_span) || pop_span == 0) pop_span <- 1
    
    for (g in unique(grp)) {
      sel <- which(grp == g & is.finite(ax))
      if (length(sel) < 5L) next
      a    <- ax[sel]
      span <- diff(range(a))
      lo   <- stats::quantile(a, q,     names = FALSE)   # extreme membership
      hi   <- stats::quantile(a, 1 - q, names = FALSE)   #   (for date pooling)
      # Normalisation anchors at the near-maximal quantiles, NOT at q/(1-q):
      # the breeding extreme is a spatial limit that may be under-represented
      # in time (skipped-breeding years sit at a lower settlement and would
      # drag a 90th-percentile anchor down, clamping the true extreme to 1).
      a_lo  <- stats::quantile(a, anchor_q,     names = FALSE)
      a_hi  <- stats::quantile(a, 1 - anchor_q, names = FALSE)
      denom <- a_hi - a_lo
      p_norm[sel] <- if (denom <= 0) 0.5 else
        pmin(1, pmax(0, (a - a_lo) / denom))
      # only well-travelled individuals vote on polarity / anchors
      if (span >= min_span_frac * pop_span) {
        is_hi <- a >= hi; is_lo <- a <= lo
        hi_dates <- c(hi_dates, doy[sel][is_hi])
        lo_dates <- c(lo_dates, doy[sel][is_lo])
        hi_lat   <- c(hi_lat,   lat[sel][is_hi])
        lo_lat   <- c(lo_lat,   lat[sel][is_lo])
      }
    }
    
    # polarity: which axis end is the breeding end?
    # If the user supplied breed_doy / winter_doy, orient the SPATIAL polarity
    # to match those dates -- the breeding end is the axis extreme whose
    # occupancy dates sit nearest the breeding date (or farthest from the
    # wintering date). This keeps the spatial and temporal terms consistent
    # even when the automatic concentration test gets the polarity wrong (the
    # override previously fixed only the temporal tent, leaving the spatial
    # term flipped -> unusable MPI). Without an override, breeding end = the
    # temporally more concentrated extreme.
    R_hi <- circ_R(hi_dates); R_lo <- circ_R(lo_dates)
    m_hi <- circ_mean(hi_dates); m_lo <- circ_mean(lo_dates)
    cdist <- function(a, b) { d <- abs((a - b) %% 365); pmin(d, 365 - d) }
    if (!is.null(breed_doy) && is.finite(m_hi) && is.finite(m_lo)) {
      breed_is_high <- cdist(m_hi, breed_doy) <= cdist(m_lo, breed_doy)
      polarity_by   <- "breed_doy override"
    } else if (!is.null(winter_doy) && is.finite(m_hi) && is.finite(m_lo)) {
      breed_is_high <- cdist(m_hi, winter_doy) >= cdist(m_lo, winter_doy)
      polarity_by   <- "winter_doy override"
    } else if (is.na(R_hi) && is.na(R_lo)) {
      breed_is_high <- TRUE                       # degenerate; arbitrary
      polarity_by   <- "degenerate"
    } else if (is.na(R_hi)) {
      breed_is_high <- FALSE; polarity_by <- "concentration"
    } else if (is.na(R_lo)) {
      breed_is_high <- TRUE;  polarity_by <- "concentration"
    } else if (abs(R_hi - R_lo) < 0.05) {
      # near-tie: fall back to the poleward (higher |lat|) end
      breed_is_high <- abs(mean(hi_lat)) >= abs(mean(lo_lat))
      polarity_by   <- "poleward (near-tie)"
    } else {
      breed_is_high <- R_hi > R_lo; polarity_by <- "concentration"
    }
    
    if (!breed_is_high) p_norm <- 1 - p_norm
    list(p = p_norm,
         breed_dates  = if (breed_is_high) hi_dates else lo_dates,
         winter_dates = if (breed_is_high) lo_dates else hi_dates,
         R_hi = R_hi, R_lo = R_lo, breed_is_high = breed_is_high,
         polarity_by = polarity_by)
  }
  
  prim <- process_axis(ax1)
  sec  <- if (use_secondary) process_axis(ax2) else NULL
  
  # ---- temporal anchors (population, two day-of-year points) -----------
  bdoy <- if (!is.null(breed_doy))  breed_doy  else circ_mean(prim$breed_dates)
  wdoy <- if (!is.null(winter_doy)) winter_doy else circ_mean(prim$winter_dates)
  if (is.na(bdoy) || is.na(wdoy))
    stop("Could not derive temporal anchors (no individual spans the ",
         "axis by `min_span_frac`). Supply `breed_doy`/`winter_doy`, or ",
         "lower `min_span_frac`.", call. = FALSE)
  
  # triangular tent d_norm over day-of-year: 0 at the wintering anchor,
  # rising linearly to 1 at the breeding anchor, then falling linearly back
  # to 0. The two arcs (spring up, autumn down) fold onto the same value, so
  # MPI is a direction-agnostic breeding-phase index.
  tent <- function(t, w, b) {
    L1 <- (b - w) %% 365; L2 <- 365 - L1
    if (L1 == 0 || L2 == 0) return(rep(0.5, length(t)))
    dist <- (t - w) %% 365
    ifelse(dist <= L1, dist / L1, 1 - (dist - L1) / L2)
  }
  
  d_norm <- tent(doy, wdoy, bdoy)
  
  # ---- combine ----------------------------------------------------------
  if (is.null(sec)) {
    mpi <- sqrt(prim$p * d_norm)
  } else {
    mpi <- (prim$p * sec$p * d_norm)^(1 / 3)
  }
  
  if (verbose) {
    message(sprintf(
      "[compute_mpi] axis: PC1 = %.1f%% var (%s)%s",
      100 * pc1_var,
      if (use_secondary) "mixed, PC1+PC2" else "unidirectional, PC1",
      ""))
    message(sprintf(
      "[compute_mpi] breed_doy = %.0f, winter_doy = %.0f",
      bdoy, wdoy))
    message(sprintf(
      "[compute_mpi] polarity: breeding = %s axis end (R_breed = %.2f vs R_winter = %.2f) [by %s]",
      if (prim$breed_is_high) "high" else "low",
      if (prim$breed_is_high) prim$R_hi else prim$R_lo,
      if (prim$breed_is_high) prim$R_lo else prim$R_hi,
      prim$polarity_by))
  }
  
  points$mpi <- mpi
  attr(points, "mpi_winter_doy") <- wdoy
  attr(points, "mpi_breed_doy")  <- bdoy
  attr(points, "mpi_pc1_var")    <- pc1_var
  points
}