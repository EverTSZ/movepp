
#' Compute step speed from successive GPS fixes
#'
#' For each point, computes a step-based speed using the great-circle (or
#' Euclidean) distance to the adjacent fix(es) in time, divided by the
#' time interval. Operates per individual when an individual identifier
#' is supplied.
#'
#' @section Why use this instead of tracker-reported speed:
#' Tracker-reported speed varies between devices: some report Doppler
#' instantaneous speed, others report displacement-over-interval, others
#' use proprietary smoothing. Many low-cost trackers do not report
#' speed at all. Computing step speed from raw fix coordinates gives a
#' reproducible, device-agnostic variable directly comparable across
#' datasets and consistent with the spatial scale of the analysis.
#'
#' @section Direction modes:
#' The `direction` argument controls how each fix's speed is defined:
#' \itemize{
#'   \item `"centered"` (default): the total distance to both adjacent
#'     fixes divided by the total time, giving a window-averaged speed
#'     around the focal fix. Symmetrically captures both take-off and
#'     landing events, and produces no edge NAs (the first and last fix
#'     degenerate to forward and backward respectively).
#'   \item `"backward"`: distance from the previous fix to the focal
#'     fix, divided by their time interval. Standard in movement
#'     ecology (e.g. \pkg{adehabitatLT}, \pkg{momentuHMM}) but biases
#'     detection toward landing events (the first fix of each
#'     individual is NA).
#'   \item `"forward"`: distance from the focal fix to the next fix,
#'     divided by their time interval. Biases detection toward take-off
#'     events (the last fix of each individual is NA).
#' }
#'
#' @param points An `sf` POINT object.
#' @param time_col Column name (string) for timestamp; must be POSIXct
#'   or Date.
#' @param individual_col Column name (string) for individual identifier,
#'   or `NULL` to treat all points as one trajectory.
#' @param unit Output unit for `step_speed`: `"km/h"` (default),
#'   `"m/s"`, or `"km/day"`.
#' @param direction Direction of the step calculation:
#'   `"centered"` (default), `"backward"`, or `"forward"`. See
#'   *Direction modes*.
#' @param max_gap_hours Numeric; if the time gap to the relevant
#'   adjacent fix exceeds this threshold, that side is treated as
#'   missing (long gaps usually indicate tracker dropouts, not actual
#'   movement). Default `Inf` (no filtering). When `direction =
#'   "centered"` the filter is applied to each side independently;
#'   the centered estimate then degenerates to whichever side remains
#'   valid (or NA if both sides are filtered).
#'
#' @return The input `sf` object with three new columns:
#'   \itemize{
#'     \item `step_distance_km`: distance used in the calculation (km).
#'       For `"centered"`, the sum of distances to both adjacent fixes.
#'     \item `step_dt_hours`: time interval used in the calculation
#'       (hours). For `"centered"`, the sum of intervals to both
#'       adjacent fixes.
#'     \item `step_speed`: `step_distance_km / step_dt_hours` in the
#'       chosen unit.
#'   }
#'
#' @examples
#' \dontrun{
#' data(godwit_demo)
#' godwit2 <- compute_step_speed(godwit_demo,
#'                               time_col       = "time",
#'                               individual_col = "individual",
#'                               unit           = "km/h",
#'                               direction      = "centered")
#' summary(godwit2$step_speed)
#'
#' # Feed into Local Moran's I
#' seg <- balm_segmentation(godwit2,
#'                           variable_col   = "step_speed",
#'                           individual_col = "individual")
#' }
#'
#' @seealso [balm_segmentation()] for using the computed step speed
#'   in spatial autocorrelation analysis.
#'
#' @export
compute_step_speed <- function(points,
                               time_col       = "time",
                               individual_col = NULL,
                               unit           = c("km/h", "m/s", "km/day"),
                               direction      = c("centered",
                                                  "backward",
                                                  "forward"),
                               max_gap_hours  = Inf) {

  if (!inherits(points, "sf"))
    stop("`points` must be an `sf` object.", call. = FALSE)
  if (!time_col %in% names(points))
    stop("Column `", time_col, "` not found in `points`.", call. = FALSE)
  if (!inherits(points[[time_col]], c("POSIXct", "POSIXt", "Date")))
    stop("`time_col` must reference a POSIXct or Date column.",
         call. = FALSE)
  unit      <- match.arg(unit)
  direction <- match.arg(direction)

  n <- nrow(points)
  step_distance_km <- rep(NA_real_, n)
  step_dt_hours    <- rep(NA_real_, n)

  times <- points[[time_col]]
  geom  <- sf::st_geometry(points)

  if (is.null(individual_col)) {
    groups <- list(seq_len(n))
  } else {
    if (!individual_col %in% names(points))
      stop("Column `", individual_col, "` not found in `points`.",
           call. = FALSE)
    groups <- split(seq_len(n), points[[individual_col]])
  }

  for (idx in groups) {
    if (length(idx) < 2L) next

    sub_order <- order(times[idx])
    sub_idx   <- idx[sub_order]

    sub_geom  <- geom[sub_idx]
    sub_times <- times[sub_idx]
    n_sub     <- length(sub_idx)

    dists_km <- as.numeric(sf::st_distance(
      sub_geom[2:n_sub],
      sub_geom[1:(n_sub - 1L)],
      by_element = TRUE
    )) / 1000
    dt_hr <- as.numeric(diff(sub_times), units = "secs") / 3600

    bwd_d <- rep(NA_real_, n_sub)
    bwd_t <- rep(NA_real_, n_sub)
    fwd_d <- rep(NA_real_, n_sub)
    fwd_t <- rep(NA_real_, n_sub)

    bwd_d[2:n_sub]        <- dists_km
    bwd_t[2:n_sub]        <- dt_hr
    fwd_d[1:(n_sub - 1L)] <- dists_km
    fwd_t[1:(n_sub - 1L)] <- dt_hr

    bwd_invalid <- !is.na(bwd_t) & bwd_t > max_gap_hours
    fwd_invalid <- !is.na(fwd_t) & fwd_t > max_gap_hours
    bwd_d[bwd_invalid] <- NA_real_
    bwd_t[bwd_invalid] <- NA_real_
    fwd_d[fwd_invalid] <- NA_real_
    fwd_t[fwd_invalid] <- NA_real_

    if (direction == "backward") {
      step_d <- bwd_d
      step_t <- bwd_t
    } else if (direction == "forward") {
      step_d <- fwd_d
      step_t <- fwd_t
    } else {
      both_na <- is.na(bwd_d) & is.na(fwd_d)
      step_d  <- ifelse(is.na(bwd_d), 0, bwd_d) +
                 ifelse(is.na(fwd_d), 0, fwd_d)
      step_t  <- ifelse(is.na(bwd_t), 0, bwd_t) +
                 ifelse(is.na(fwd_t), 0, fwd_t)
      step_d[both_na] <- NA_real_
      step_t[both_na] <- NA_real_
    }

    zero_t <- !is.na(step_t) & step_t == 0
    step_t[zero_t] <- NA_real_

    step_distance_km[sub_idx] <- step_d
    step_dt_hours[sub_idx]    <- step_t
  }

  step_speed_kmh <- step_distance_km / step_dt_hours
  step_speed <- switch(unit,
    "km/h"   = step_speed_kmh,
    "m/s"    = step_speed_kmh * 1000 / 3600,
    "km/day" = step_speed_kmh * 24
  )

  points$step_distance_km <- step_distance_km
  points$step_dt_hours    <- step_dt_hours
  points$step_speed       <- step_speed

  points
}

