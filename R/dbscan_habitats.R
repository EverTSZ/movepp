#' Delineate temporary habitats by DBSCAN at a behaviourally-derived scale
#'
#' Clusters stationary tracking fixes into temporary habitats using standard
#' DBSCAN (Ester et al. 1996), applied independently to each individual. Unlike
#' density-hierarchy methods (HDBSCAN), both parameters here are fixed and
#' interpretable rather than auto-extracted:
#'
#' * `eps` is a fixed, biologically meaningful spatial scale -- the
#'   species-typical maximum within-habitat movement, obtained from the
#'   population step-length distribution (the upper mode of local, non-flight
#'   steps). Occupancy at this scale is self-similar (scale-free), so a single
#'   explicit scale is more defensible than a brittle auto-selected optimum.
#' * `minPts` is a unified noise / minimum-duration floor -- the knee of the
#'   cluster-count vs `minPts` curve, which for the Black-tailed Godwit data
#'   coincides with a minimum residence of ~3 days at a 6-hour sampling
#'   interval (`minPts = 12`).
#'
#' Time is deliberately excluded: habitats are defined purely from spatial
#' occupancy, leaving temporal use to be analysed downstream.
#'
#' @param x An sf POINT object, typically the stationary (LL/LH) subset from
#'   [balm_segmentation()].
#' @param individual_col Name of the column identifying individuals. Clustering
#'   is performed independently within each individual. If `NULL`, all rows are
#'   treated as a single individual.
#' @param eps Neighbourhood radius in kilometres. **Required, no default** --
#'   it must be derived from the data, typically via [detect_habitat_params()],
#'   not hard-coded. Two fixes join the same habitat when reachable through hops
#'   no larger than `eps`.
#' @param minPts Minimum number of points within `eps` for a core point.
#'   **Required, no default** -- typically from [detect_habitat_params()]. Acts
#'   as a uniform density / noise floor; interpretable as a minimum residence
#'   time given the sampling interval.
#' @param drop_noise Logical; if `TRUE` (default) noise points (DBSCAN cluster
#'   `0`) are dropped from the returned object.
#'
#' @return The input sf object with an added integer column `cluster_id`
#'   giving the per-individual habitat id (`1..k`; `0` = noise, retained only
#'   when `drop_noise = FALSE`).
#'
#' @references Ester, M., Kriegel, H.-P., Sander, J., & Xu, X. (1996).
#'   A density-based algorithm for discovering clusters in large spatial
#'   databases with noise. *Proc. KDD-96*, 226-231.
#'
#' @examples
#' \dontrun{
#' stationary <- seg[seg$cluster_type %in% c("LL", "LH"), ]
#' ## 1. derive thresholds from the data (never hard-code them)
#' p <- detect_habitat_params(stationary, individual_col = "Individual",
#'                            time_col = "Time")
#' ## 2. run transparent per-individual DBSCAN with the derived thresholds
#' hab <- dbscan_habitats(stationary, individual_col = "Individual",
#'                        eps = p$eps, minPts = p$minPts)
#' tapply(hab$cluster_id, hab$Individual, function(z) length(unique(z)))
#' }
#' @seealso [detect_habitat_params()] to derive `eps` and `minPts` from data.
#' @export
dbscan_habitats <- function(x, individual_col = NULL, eps, minPts,
                            drop_noise = TRUE) {
  if (!inherits(x, "sf"))
    stop("`x` must be an sf object.", call. = FALSE)
  if (missing(eps) || missing(minPts))
    stop("`eps` and `minPts` must be supplied (e.g. from ",
         "detect_habitat_params()); they are intentionally not defaulted.",
         call. = FALSE)
  if (!requireNamespace("dbscan", quietly = TRUE))
    stop("Package 'dbscan' is required for dbscan_habitats().", call. = FALSE)
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for dbscan_habitats().", call. = FALSE)

  minPts <- as.integer(minPts)
  grp <- if (is.null(individual_col)) rep(1L, nrow(x)) else x[[individual_col]]
  if (is.null(individual_col) && nrow(x) == 0L)
    return(x)

  cluster_id <- integer(nrow(x))   # 0 = noise

  for (g in unique(grp)) {
    sel <- which(grp == g)
    if (length(sel) < minPts) next                      # too few fixes to seed a habitat
    co   <- sf::st_coordinates(x[sel, ])
    latm <- mean(co[, 2]); rad <- pi / 180
    # local equirectangular projection to kilometres (cos-latitude scaling)
    xy   <- cbind(co[, 1] * cos(latm * rad) * 111.320,
                  co[, 2] * 111.320)
    cluster_id[sel] <- dbscan::dbscan(xy, eps = eps, minPts = minPts)$cluster
  }

  x$cluster_id <- cluster_id
  if (drop_noise) x <- x[x$cluster_id > 0, ]
  x
}
