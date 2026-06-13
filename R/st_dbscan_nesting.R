#' Detect Bird Nesting Events with Stratified Spatiotemporal DBSCAN
#'
#' Identifies nesting locations and hatching periods from high-resolution
#' tracking data using a three-tier stratified ST-DBSCAN scheme. The
#' detection rests on the assumption that a true nest yields GPS fixes
#' approximately bivariate-normally distributed around the actual nest
#' location: 95% of points within \eqn{2\sigma}, 68% within \eqn{1\sigma},
#' and ~38% within \eqn{0.5\sigma} of the centre. A point is classified
#' as a nest point only if it survives the density threshold at all
#' three sigma tiers simultaneously, yielding markedly fewer false
#' positives than single-tier DBSCAN.
#'
#' @section Algorithm:
#' For each individual (or pooled dataset if `individual_col = NULL`):
#' \enumerate{
#'   \item Run ST-DBSCAN at three nested density tiers:
#'     \itemize{
#'       \item \strong{2\eqn{\sigma}}: `min_pts` x 1.0,
#'         `search_dist` x 1.0.
#'       \item \strong{1\eqn{\sigma}}: `min_pts` x 0.68,
#'         `search_dist` x 0.5.
#'       \item \strong{0.5\eqn{\sigma}}: `min_pts` x 0.38,
#'         `search_dist` x 0.25.
#'     }
#'   \item Retain only points classified into a non-noise cluster in
#'     \strong{all three tiers} (membership gate against false positives).
#'   \item Assign each surviving point the cluster label of the
#'     \strong{loosest (2\eqn{\sigma}) tier}, whose large radius collapses
#'     one nest's entire fix cloud into a single `nest_id`. (Labelling by the
#'     strictest tier would split one real nest into several sub-clusters,
#'     producing many spurious nest records for a single nest.)
#'   \item Compute the mean centre, first/last timestamp, and mean
#'     timestamp of each detected nest.
#' }
#'
#' @section Coordinate units:
#' `search_dist` is interpreted in the linear units of the input CRS.
#' If `points` is in WGS84 (EPSG:4326), distances are in degrees and
#' the function will warn. For accurate nest detection at scales of
#' tens of metres, project to a metric CRS first (e.g., a local UTM
#' zone) and supply `search_dist` in metres.
#'
#' @param points An `sf` POINT object with high-resolution GPS fixes
#'   during the breeding season.
#' @param time_col Column name (string) of a POSIXct/Date time column.
#' @param individual_col Column name (string) of the individual ID
#'   (default `NULL` = pool all points).
#' @param min_pts Integer; baseline minimum points per cluster at the
#'   2\eqn{\sigma} tier (default 40).
#' @param search_dist Numeric; baseline spatial radius at the
#'   2\eqn{\sigma} tier, in CRS linear units (default 41, typical for
#'   2 x 20m GPS error).
#' @param time_window Character or `difftime`; the temporal window
#'   used by ST-DBSCAN (default `"3 weeks"`, suitable for a single
#'   nesting attempt).
#' @param sigma_ratios Length-3 numeric vector of (min_pts, dist)
#'   multipliers for the (2\eqn{\sigma}, 1\eqn{\sigma}, 0.5\eqn{\sigma})
#'   tiers. Default `c(1.0, 0.68, 0.38)` for min_pts and is paired
#'   internally with `c(1.0, 0.5, 0.25)` for distance.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with two `sf` elements:
#'   \itemize{
#'     \item `points`: the input `sf`, augmented with `nest_id`
#'       (`NA` if the point is not a nest point) and `near_dist`
#'       (distance from the point to its nest centre).
#'     \item `nests`: one row per detected nest, with columns
#'       `nest_id`, `individual`, `n_points`, `start_time`,
#'       `end_time`, `mean_time`, and geometry (the mean centre).
#'   }
#'
#' @references
#' Birant, D. & Kut, A. (2007). ST-DBSCAN: An algorithm for clustering
#' spatial-temporal data. *Data & Knowledge Engineering* 60(1): 208-221.
#' \doi{10.1016/j.datak.2006.01.013}
#'
#' @examples
#' \dontrun{
#' # Pied Avocet nesting data (project to UTM first for metric distance)
#' nests <- detect_nests(avocet_breeding_pts,
#'                       time_col = "time",
#'                       individual_col = "individual",
#'                       min_pts = 40,
#'                       search_dist = 41,
#'                       time_window = "3 weeks")
#' nests$nests  # mean centre + time bounds of each nest
#' }
#'
#' @export
detect_nests <- function(points,
                         time_col,
                         individual_col = NULL,
                         min_pts        = 40L,
                         search_dist    = 41,
                         time_window    = "3 weeks",
                         sigma_ratios   = c(1.0, 0.68, 0.38),
                         verbose        = TRUE) {
  # ---- Input validation -------------------------------------------------
  if (!inherits(points, "sf")) {
    stop("`points` must be an `sf` object.", call. = FALSE)
  }
  if (!time_col %in% names(points)) {
    stop("Column `", time_col, "` not found.", call. = FALSE)
  }
  if (!inherits(points[[time_col]], c("POSIXct", "POSIXt", "Date"))) {
    stop("`time_col` must be a POSIXct or Date column.", call. = FALSE)
  }
  if (length(sigma_ratios) != 3L) {
    stop("`sigma_ratios` must be a length-3 numeric vector.", call. = FALSE)
  }
  # Warn about geographic CRS
  crs_info <- sf::st_crs(points)
  if (!is.na(crs_info) && isTRUE(crs_info$IsGeographic)) {
    warning("Input is in a geographic CRS; `search_dist` will be ",
            "interpreted as degrees. Project to a metric CRS for ",
            "accurate nest detection.", call. = FALSE)
  }
  # Parse time window
  time_window_secs <- .parse_time_window(time_window)
  # ---- Setup tier parameters --------------------------------------------
  dist_ratios <- c(1.0, 0.5, 0.25)
  tier_params <- data.frame(
    name    = c("2sigma", "1sigma", "0.5sigma"),
    min_pts = as.integer(round(min_pts * sigma_ratios)),
    dist    = search_dist * dist_ratios
  )
  if (verbose) {
    message("[detect_nests] Tier parameters:")
    for (i in seq_len(nrow(tier_params))) {
      message(sprintf("  %s: min_pts=%d, dist=%.2f",
                      tier_params$name[i],
                      tier_params$min_pts[i],
                      tier_params$dist[i]))
    }
  }
  # ---- Group iteration --------------------------------------------------
  if (is.null(individual_col)) {
    groups <- list("all" = seq_len(nrow(points)))
  } else {
    groups <- split(seq_len(nrow(points)),
                    as.character(points[[individual_col]]))
  }
  points$nest_id   <- rep(NA_integer_, nrow(points))
  points$near_dist <- rep(NA_real_, nrow(points))
  all_nests <- list()
  next_nest_id <- 1L
  for (g_name in names(groups)) {
    g_idx <- groups[[g_name]]
    n_pts <- length(g_idx)
    if (n_pts < tier_params$min_pts[1]) {
      if (verbose) {
        message(sprintf("[detect_nests] Skip '%s' (%d < min_pts)",
                        g_name, n_pts))
      }
      next
    }
    if (verbose) {
      message(sprintf("[detect_nests] Processing '%s' (%d points)...",
                      g_name, n_pts))
    }
    g_pts    <- points[g_idx, ]
    g_coords <- sf::st_coordinates(g_pts)
    g_times  <- g_pts[[time_col]]
    # Run ST-DBSCAN at each tier
    tier_clusters <- vector("list", 3L)
    for (t in 1:3) {
      tier_clusters[[t]] <- .st_dbscan_core(
        coords        = g_coords,
        times         = g_times,
        eps_spatial   = tier_params$dist[t],
        eps_temporal  = time_window_secs,
        min_pts       = tier_params$min_pts[t]
      )$cluster
    }
    # Membership gate: a point counts as a nest point only if it is in a
    # non-noise cluster in ALL three tiers (dense at every scale, i.e. a
    # bivariate-normal peak rather than a flat aggregation).
    in_all_tiers <- (tier_clusters[[1]] > 0L) &
      (tier_clusters[[2]] > 0L) &
      (tier_clusters[[3]] > 0L)
    if (!any(in_all_tiers)) {
      if (verbose) {
        message(sprintf("[detect_nests] No nest detected for '%s'.",
                        g_name))
      }
      next
    }
    # Identity: label the survivors by the LOOSEST tier (2 sigma). Its large
    # radius collapses one nest's whole fix cloud into a single cluster, so
    # one real nest yields exactly one nest_id. (The 0.5 sigma tier would
    # fragment a single nest into several sub-clusters -> many spurious
    # records.) The three-tier intersection above is only the gate.
    nest_labels <- tier_clusters[[1]]
    nest_labels[!in_all_tiers] <- 0L
    # Re-number to globally unique nest IDs
    raw_nest_ids <- unique(nest_labels[nest_labels > 0L])
    id_map <- stats::setNames(
      seq_along(raw_nest_ids) + (next_nest_id - 1L),
      as.character(raw_nest_ids)
    )
    nest_ids <- rep(NA_integer_, n_pts)
    valid_idx <- nest_labels > 0L
    nest_ids[valid_idx] <-
      as.integer(id_map[as.character(nest_labels[valid_idx])])
    points$nest_id[g_idx] <- nest_ids
    next_nest_id <- next_nest_id + length(raw_nest_ids)
    # Compute nest centres for this group
    for (nid in unique(nest_ids[!is.na(nest_ids)])) {
      member_idx <- which(nest_ids == nid)
      member_coords <- g_coords[member_idx, , drop = FALSE]
      member_times  <- g_times[member_idx]
      mean_lon <- mean(member_coords[, "X"])
      mean_lat <- mean(member_coords[, "Y"])
      # Distance from each member to centre
      dists <- sqrt((member_coords[, "X"] - mean_lon)^2 +
                      (member_coords[, "Y"] - mean_lat)^2)
      points$near_dist[g_idx[member_idx]] <- dists
      nest_row <- data.frame(
        nest_id    = nid,
        individual = g_name,
        n_points   = length(member_idx),
        start_time = min(member_times),
        end_time   = max(member_times),
        mean_time  = mean(member_times),
        lon        = mean_lon,
        lat        = mean_lat
      )
      all_nests[[length(all_nests) + 1L]] <- nest_row
    }
  }
  # ---- Assemble nest output ---------------------------------------------
  if (length(all_nests) == 0L) {
    if (verbose) message("[detect_nests] No nests detected overall.")
    nests_sf <- sf::st_sf(
      nest_id    = integer(0),
      individual = character(0),
      n_points   = integer(0),
      start_time = as.POSIXct(character(0)),
      end_time   = as.POSIXct(character(0)),
      mean_time  = as.POSIXct(character(0)),
      geometry   = sf::st_sfc(crs = crs_info)
    )
  } else {
    nests_df <- do.call(rbind, all_nests)
    nests_sf <- sf::st_as_sf(nests_df, coords = c("lon", "lat"),
                             crs = crs_info, remove = FALSE)
  }
  if (verbose) {
    n_nest <- if (inherits(nests_sf, "sf")) nrow(nests_sf) else 0L
    message(sprintf("[detect_nests] Done: %d nest(s) detected across %d group(s).",
                    n_nest, length(groups)))
  }
  list(points = points, nests = nests_sf)
}

#' Internal: Core ST-DBSCAN clustering algorithm
#'
#' Implements the spatiotemporal DBSCAN of Birant & Kut (2007).
#' Not exported; called by [detect_nests()].
#'
#' @param coords Numeric matrix or 2-column data frame of (x, y).
#' @param times POSIXct/Date vector aligned to `coords` rows.
#' @param eps_spatial Numeric; spatial neighbourhood radius (CRS units).
#' @param eps_temporal Numeric; temporal neighbourhood radius (seconds).
#' @param min_pts Integer; minimum points per cluster.
#'
#' @return A list with `cluster` (integer, 0 = noise, 1..k = clusters)
#'   and `n_clusters`.
#'
#' @keywords internal
.st_dbscan_core <- function(coords, times,
                            eps_spatial, eps_temporal, min_pts) {
  n <- nrow(coords)
  cluster <- integer(n)   # 0 = noise/unassigned
  visited <- logical(n)
  cid     <- 0L
  t_sec <- as.numeric(times)
  for (i in seq_len(n)) {
    if (visited[i]) next
    visited[i] <- TRUE
    # Find neighbours of i (spatial AND temporal within bounds)
    d_spatial <- sqrt((coords[, 1] - coords[i, 1])^2 +
                        (coords[, 2] - coords[i, 2])^2)
    d_temporal <- abs(t_sec - t_sec[i])
    nbrs <- which(d_spatial <= eps_spatial &
                    d_temporal <= eps_temporal)
    if (length(nbrs) < min_pts) {
      next  # noise (remains 0)
    }
    cid <- cid + 1L
    cluster[i] <- cid
    # BFS expansion
    seed <- setdiff(nbrs, i)
    while (length(seed) > 0L) {
      j <- seed[1L]
      seed <- seed[-1L]
      if (!visited[j]) {
        visited[j] <- TRUE
        dj_s <- sqrt((coords[, 1] - coords[j, 1])^2 +
                       (coords[, 2] - coords[j, 2])^2)
        dj_t <- abs(t_sec - t_sec[j])
        nbrs_j <- which(dj_s <= eps_spatial & dj_t <= eps_temporal)
        if (length(nbrs_j) >= min_pts) {
          seed <- union(seed, nbrs_j)
        }
      }
      if (cluster[j] == 0L) {
        cluster[j] <- cid
      }
    }
  }
  list(cluster = cluster, n_clusters = cid)
}

#' Internal: Parse a human-readable time window to seconds
#'
#' @param tw Character (e.g., `"3 weeks"`, `"1 day"`) or a `difftime`.
#'
#' @return Numeric, number of seconds.
#'
#' @keywords internal
.parse_time_window <- function(tw) {
  if (inherits(tw, "difftime")) {
    return(as.numeric(tw, units = "secs"))
  }
  if (!is.character(tw)) {
    stop("`time_window` must be a character string or difftime.",
         call. = FALSE)
  }
  parts <- strsplit(trimws(tw), "\\s+")[[1]]
  if (length(parts) != 2L) {
    stop("`time_window` must look like '3 weeks' or '1 day'.",
         call. = FALSE)
  }
  num  <- suppressWarnings(as.numeric(parts[1]))
  unit <- tolower(parts[2])
  if (is.na(num)) {
    stop("Cannot parse number from '", tw, "'.", call. = FALSE)
  }
  mult <- switch(
    sub("s$", "", unit),
    second = 1, sec = 1,
    minute = 60, min = 60,
    hour   = 3600,
    day    = 86400,
    week   = 86400 * 7,
    month  = 86400 * 30,
    year   = 86400 * 365,
    stop("Unknown time unit '", unit, "'.", call. = FALSE)
  )
  num * mult
}