#' Compute Solar Elevation Angle at Each Tracking Point
#'
#' Calculates the solar elevation angle (in degrees, with negative
#' values below the horizon) for each tracking point, given its
#' coordinates and timestamp. Used downstream by [classify_circadian()]
#' to partition activity into day, twilight, and night periods.
#'
#' @param points An `sf` POINT object.
#' @param time_col Column name (string) of a POSIXct/POSIXt time
#'   column. Must include time zone information for correct local-time
#'   computation.
#' @param verbose Logical; print progress.
#'
#' @return The input `sf` with a new `solar_elevation` numeric column
#'   in degrees (positive = above horizon, negative = below).
#'
#' @references
#' Whitener, B. (2025). *Pysolar: Python Library for Solar Position
#' Calculations*. The R equivalent computation is performed by the
#' `suncalc` package (Thieurmel & Elmarhraoui).
#'
#' @examples
#' \dontrun{
#' demo <- make_demo_track(n_individuals = 1, n_points = 100)
#' demo <- compute_solar_elevation(demo, time_col = "time")
#' hist(demo$solar_elevation)
#' }
#'
#' @export
compute_solar_elevation <- function(points,
                                    time_col,
                                    verbose = TRUE) {

  if (!inherits(points, "sf")) {
    stop("`points` must be an `sf` object.", call. = FALSE)
  }
  if (!time_col %in% names(points)) {
    stop("Column `", time_col, "` not found.", call. = FALSE)
  }
  if (!requireNamespace("suncalc", quietly = TRUE)) {
    stop("Package 'suncalc' is required.", call. = FALSE)
  }

  tvec <- points[[time_col]]
  if (!inherits(tvec, c("POSIXct", "POSIXt"))) {
    stop("`time_col` must reference a POSIXct/POSIXt column.",
         call. = FALSE)
  }

  # Get coords as lon/lat (suncalc expects geographic CRS)
  cur_crs <- sf::st_crs(points)
  if (is.na(cur_crs) || !isTRUE(cur_crs$IsGeographic)) {
    if (verbose) {
      message("[compute_solar_elevation] Transforming coordinates to ",
              "EPSG:4326 for solar position computation.")
    }
    pts_geo <- sf::st_transform(points, 4326)
  } else {
    pts_geo <- points
  }
  coords <- sf::st_coordinates(pts_geo)

  if (verbose) {
    message(sprintf("[compute_solar_elevation] Computing for %d points...",
                    nrow(points)))
  }

  # suncalc::getSunlightPosition returns altitude in radians
  # suncalc requires a data.frame when lat/lon vary per row
  sun_pos <- suncalc::getSunlightPosition(
    data = data.frame(
      date = tvec,
      lat  = coords[, "Y"],
      lon  = coords[, "X"]
    )
  )

  # Convert to degrees
  points$solar_elevation <- sun_pos$altitude * 180 / pi

  if (verbose) {
    message(sprintf(
      "[compute_solar_elevation] Done. Range: %.1f to %.1f deg",
      min(points$solar_elevation, na.rm = TRUE),
      max(points$solar_elevation, na.rm = TRUE)))
  }

  points
}


#' Classify Tracking Points into Day / Twilight / Night
#'
#' Partitions tracking points by solar elevation angle into three
#' circadian categories: **Day** (sun above horizon), **Twilight**
#' (civil twilight, sun between the two thresholds), and **Night**
#' (sun well below horizon). The default thresholds follow the standard
#' definition of civil twilight (0 deg and -6 deg).
#'
#' @param points An `sf` object with a solar elevation column (typically
#'   the output of [compute_solar_elevation()]).
#' @param elevation_col Column name (string) of solar elevation in
#'   degrees (default `"solar_elevation"`).
#' @param day_threshold Numeric; angle (deg) above which a point is
#'   classified as Day (default 0).
#' @param twilight_threshold Numeric; angle (deg) above which a point
#'   is classified as Twilight (default -6, the civil-twilight bound).
#' @param verbose Logical; print summary.
#'
#' @return The input `sf` with a new `circadian` factor column
#'   (levels: `"Day"`, `"Twilight"`, `"Night"`).
#'
#' @examples
#' \dontrun{
#' demo <- compute_solar_elevation(demo, time_col = "time")
#' demo <- classify_circadian(demo)
#' table(demo$circadian)
#' }
#'
#' @export
classify_circadian <- function(points,
                               elevation_col      = "solar_elevation",
                               day_threshold      = 0,
                               twilight_threshold = -6,
                               verbose            = TRUE) {

  if (!inherits(points, "sf")) {
    stop("`points` must be an `sf` object.", call. = FALSE)
  }
  if (!elevation_col %in% names(points)) {
    stop("Column `", elevation_col, "` not found. ",
         "Did you run `compute_solar_elevation()` first?",
         call. = FALSE)
  }
  if (twilight_threshold >= day_threshold) {
    stop("`twilight_threshold` must be < `day_threshold`.",
         call. = FALSE)
  }

  elev <- points[[elevation_col]]

  category <- ifelse(is.na(elev), NA_character_,
                ifelse(elev > day_threshold, "Day",
                ifelse(elev > twilight_threshold, "Twilight", "Night")))

  points$circadian <- factor(category,
                             levels = c("Day", "Twilight", "Night"))

  if (verbose) {
    tbl <- table(points$circadian, useNA = "ifany")
    message("[classify_circadian] Counts: ",
            paste(names(tbl), tbl, sep = "=", collapse = ", "))
  }

  points
}


#' Type each temporary habitat by its dominant diel period
#'
#' Aggregates the per-fix circadian labels from [classify_circadian()] to
#' the temporary-habitat (cluster) level. Each habitat is assigned the
#' diel period holding a majority (greater than `majority`) of its fixes,
#' or `"Mixed"` when no period does. The per-fix day/twilight/night label
#' is only an intermediate; the habitat-level functional type -- a daytime
#' foraging site, a night roost, a twilight transit area -- is the
#' ecological object this returns.
#'
#' @param points An `sf` object carrying a per-fix circadian label
#'   (typically the `circadian` column from [classify_circadian()]) and a
#'   temporary-habitat cluster id (the `cluster_id` from
#'   [dbscan_habitats()]).
#' @param circadian_col Per-fix diel-label column (default `"circadian"`).
#' @param cluster_col Temporary-habitat cluster-id column (default
#'   `"cluster_id"`).
#' @param individual_col Individual-id column, or `NULL` (default) to treat
#'   `cluster_col` as globally unique. When supplied, habitats are typed
#'   within each `(individual, cluster)` pair.
#' @param majority Numeric in \[0, 1\]; the minimum fraction the dominant
#'   period must reach for the habitat to take its label rather than
#'   `"Mixed"` (default 0.5).
#' @param periods Character vector of diel levels to tabulate and rank
#'   (default `c("Day", "Twilight", "Night")`).
#' @param verbose Logical; print the resulting habitat-type counts.
#'
#' @return The input `sf` with an added factor column `habitat_circadian`
#'   (levels `periods` plus `"Mixed"`), constant within each habitat. A
#'   one-row-per-habitat summary -- cluster id, fix count, per-period
#'   proportions, and assigned type -- is attached as the attribute
#'   `"habitat_circadian"`.
#'
#' @examples
#' \dontrun{
#' pts <- compute_solar_elevation(habitat_pts, time_col = "time")
#' pts <- classify_circadian(pts)
#' pts <- classify_habitat_circadian(pts, individual_col = "individual")
#' attr(pts, "habitat_circadian")
#' }
#'
#' @seealso [classify_circadian()] for the per-fix labels this aggregates.
#' @export
classify_habitat_circadian <- function(points,
                                       circadian_col  = "circadian",
                                       cluster_col    = "cluster_id",
                                       individual_col = NULL,
                                       majority       = 0.5,
                                       periods        = c("Day", "Twilight", "Night"),
                                       verbose        = TRUE) {

  if (!inherits(points, "sf"))
    stop("`points` must be an `sf` object.", call. = FALSE)
  for (col in c(circadian_col, cluster_col))
    if (!col %in% names(points))
      stop("Column `", col, "` not found in `points`.", call. = FALSE)
  if (!is.null(individual_col) && !individual_col %in% names(points))
    stop("Column `", individual_col, "` not found in `points`.", call. = FALSE)
  if (majority < 0 || majority > 1)
    stop("`majority` must be in [0, 1].", call. = FALSE)

  ind <- if (is.null(individual_col)) rep("", nrow(points))
         else as.character(points[[individual_col]])
  cl  <- as.character(points[[cluster_col]])
  key <- if (is.null(individual_col)) cl else paste(ind, cl, sep = "::")
  cls <- factor(as.character(points[[circadian_col]]), levels = periods)
  ok  <- !is.na(key) & !is.na(cls)

  lev <- c(periods, "Mixed")
  if (!any(ok)) {
    points$habitat_circadian <- factor(rep(NA_character_, nrow(points)),
                                       levels = lev)
    attr(points, "habitat_circadian") <- data.frame()
    if (verbose) message("[classify_habitat_circadian] no valid fixes to type.")
    return(points)
  }

  tab   <- table(key[ok], cls[ok])             # habitats x periods
  prop  <- prop.table(tab, 1)                  # per-habitat diel proportions
  htype <- apply(prop, 1, function(p)
    if (max(p) > majority) names(p)[which.max(p)] else "Mixed")

  points$habitat_circadian <- factor(unname(htype[key]), levels = lev)

  keys  <- rownames(prop)
  first <- match(keys, key)                    # recover ids without splitting
  summary_df <- data.frame(cluster_id = cl[first],
                           n           = as.integer(rowSums(tab)),
                           stringsAsFactors = FALSE)
  if (!is.null(individual_col))
    summary_df <- cbind(individual = ind[first], summary_df,
                        stringsAsFactors = FALSE)
  for (pd in periods)
    summary_df[[paste0("prop_", pd)]] <- round(as.numeric(prop[, pd]), 3)
  summary_df$type <- unname(htype[keys])
  rownames(summary_df) <- NULL
  attr(points, "habitat_circadian") <- summary_df

  if (verbose) {
    tt <- table(points$habitat_circadian, useNA = "ifany")
    message("[classify_habitat_circadian] habitat types: ",
            paste(names(tt), tt, sep = "=", collapse = ", "))
  }
  points
}


