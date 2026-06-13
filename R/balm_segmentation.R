#' Behaviorally-Anchored Local Moran (BALM) movement-state segmentation
#'
#' Classifies each tracking fix into a movement state -- migration,
#' stationary, or transitional -- with a behaviorally-anchored variant of
#' the Local Moran scatterplot, using only point locations and a
#' movement-rate mark (no direction, no raw temporal sequence).
#'
#' @section Method:
#' BALM adapts the Moran scatterplot (own value vs spatial lag) in two
#' ways suited to movement data:
#' \enumerate{
#'   \item \strong{Behavioral anchoring.} The high/low reference is not the
#'     arithmetic mean (standard Local Moran centering) but a data-driven
#'     behavioral threshold `c`: the onset of the flight mode. It is
#'     estimated once from the pooled (log) movement-rate of ALL individuals
#'     by Gaussian-mixture decomposition (component count chosen by BIC);
#'     `c` is the lower 5\% quantile of the fastest-mean component. Pooling
#'     yields a stable species-level threshold, applied to every individual,
#'     while the spatial neighbourhood stays per individual.
#'   \item \strong{Deterministic classification.} Each fix is labelled by
#'     the signs of its deviation `d_i = x_i - c` and its spatial lag
#'     `Sum_j w_ij d_j`; no permutation significance is applied. Spatial
#'     support is arbitrated downstream by density clustering
#'     ([dbscan_habitats()]).
#' }
#' The four states follow the Moran-scatterplot quadrants:
#' \itemize{
#'   \item `HH`: fast fix among fast neighbours -> migration.
#'   \item `LL`: slow fix among slow neighbours -> stationary site.
#'   \item `HL`: fast fix among slow neighbours -> local movement within a
#'     site (e.g. commuting or foraging burst).
#'   \item `LH`: slow fix among fast neighbours -> brief touch-down in the
#'     migratory corridor (transient stopover).
#' }
#' Because it never uses direction or fix ordering, BALM is invariant to
#' the temporal sequence and robust to sampling-rate heterogeneity.
#'
#' @param points An `sf` POINT object (geographic CRS recommended).
#' @param variable_col Column name (string) of the movement-rate
#'   attribute, e.g. `"step_speed"` from [compute_step_speed()].
#' @param individual_col Column name (string) of the individual ID, or
#'   `NULL` (default) to treat all points as one group. The behavioral
#'   reference is pooled across all individuals; the neighbourhood is
#'   computed per group. Pool one species at a time.
#' @param reference Either `"auto"` (default; species-level flight onset
#'   estimated from all individuals pooled) or a single numeric used as a
#'   fixed reference for all groups.
#' @param k Integer; nearest neighbours for the spatial weights
#'   (default 8; inverse-distance, row-standardised).
#' @param kde_adjust Numeric; passed to the flight-onset estimator.
#' @param verbose Logical; print the pooled reference and per-group counts.
#'
#' @return The input `sf` with `balm_deviation`, `balm_lag`, and a
#'   `cluster_type` factor (levels `HH`, `LL`, `HL`, `LH`; NA where
#'   `variable_col` is missing). The reference is attached as attribute
#'   `"balm_reference"`.
#'
#' @seealso [compute_step_speed()] for the input rate;
#'   [dbscan_habitats()] for downstream habitat delineation.
#'
#' @export
balm_segmentation <- function(points,
                              variable_col,
                              individual_col = NULL,
                              reference = "auto",
                              k = 8L,
                              kde_adjust = 1,
                              verbose = TRUE) {
  if (!inherits(points, "sf"))
    stop("`points` must be an `sf` object.", call. = FALSE)
  if (!variable_col %in% names(points))
    stop("Column `", variable_col, "` not found in `points`.", call. = FALSE)
  if (!is.null(individual_col) && !individual_col %in% names(points))
    stop("Column `", individual_col, "` not found in `points`.", call. = FALSE)
  if (!requireNamespace("dbscan", quietly = TRUE))
    stop("Package 'dbscan' is required.", call. = FALSE)
  k <- as.integer(k)
  n <- nrow(points)
  points$balm_deviation <- NA_real_
  points$balm_lag        <- NA_real_
  ct <- rep(NA_character_, n)
  
  if (is.null(individual_col)) {
    groups <- list(all = seq_len(n))
  } else {
    groups <- split(seq_len(n), points[[individual_col]])
  }
  refs  <- stats::setNames(rep(NA_real_, length(groups)), names(groups))
  is_ll <- isTRUE(sf::st_is_longlat(points))
  
  # ---- One species-level flight onset, pooled across ALL individuals ----
  # The threshold is a species behavioral constant, so pooling stabilises it.
  # Only the reference is pooled; the spatial neighbourhood (kNN) stays per group.
  if (is.character(reference)) {
    xall <- points[[variable_col]]
    cc <- .balm_flight_onset(xall[is.finite(xall)], adjust = kde_adjust)
    if (is.na(cc))
      stop("Flight-onset could not be detected from the pooled rates.",
           call. = FALSE)
  } else {
    cc <- as.numeric(reference)
  }
  if (verbose) message(sprintf("[balm] pooled reference c = %.3f", cc))
  
  for (gi in seq_along(groups)) {
    idx <- groups[[gi]]; gname <- names(groups)[gi]
    x <- points[[variable_col]][idx]; ok <- is.finite(x)
    if (sum(ok) <= k + 1L) {
      if (verbose) message(sprintf("[balm] '%s': too few valid points (%d); skipped.",
                                   gname, sum(ok)))
      next
    }
    idv <- idx[ok]; xv <- x[ok]
    refs[gi] <- cc
    
    sub <- points[idv, ]
    co  <- sf::st_coordinates(sub)
    if (is_ll) co[, 1] <- co[, 1] * cos(mean(co[, 2]) * pi / 180)
    knn <- dbscan::kNN(co, k = k)
    w   <- 1 / (knn$dist + 1e-9); w <- w / rowSums(w)
    
    d   <- xv - cc
    lag <- rowSums(w * matrix(d[knn$id], nrow = length(idv)))
    quad <- ifelse(d > 0 & lag > 0, "HH",
                   ifelse(d < 0 & lag < 0, "LL",
                          ifelse(d > 0 & lag < 0, "HL", "LH")))
    
    points$balm_deviation[idv] <- d
    points$balm_lag[idv]       <- lag
    ct[idv] <- quad
    
    if (verbose) {
      tb <- table(quad)
      message(sprintf("[balm] '%s': %s", gname,
                      paste(names(tb), tb, sep = "=", collapse = " ")))
    }
  }
  points$cluster_type <- factor(ct, levels = c("HH", "LL", "HL", "LH"))
  attr(points, "balm_reference") <- refs
  points
}

#' @keywords internal
.balm_flight_onset <- function(v, adjust = 1, lower_q = 0.05, G = 1:9) {
  lv <- log(v[is.finite(v) & v > 0])
  if (length(lv) < 10L) return(NA_real_)
  if (!requireNamespace("mclust", quietly = TRUE))
    stop("Package 'mclust' is required.", call. = FALSE)
  
  # Reproducible GMM (restore the caller's RNG state on exit).
  ss <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(ss)) assign(".Random.seed", ss, envir = .GlobalEnv), add = TRUE)
  set.seed(1L)
  
  fit <- tryCatch(mclust::Mclust(sort(lv), G = G, verbose = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(exp(stats::quantile(lv, lower_q, names = FALSE)))
  
  mu  <- fit$parameters$mean
  vv  <- fit$parameters$variance
  sds <- if (!is.null(vv$sigmasq)) {
    if (length(vv$sigmasq) == 1L) rep(sqrt(vv$sigmasq), length(mu)) else sqrt(vv$sigmasq)
  } else sqrt(as.numeric(vv$sigma))
  
  j <- which.max(mu)                                     # rightmost (fastest) = flight
  exp(stats::qnorm(lower_q, mean = mu[j], sd = sds[j]))  # its 5% quantile = flight onset
}