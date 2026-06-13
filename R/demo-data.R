#' Generate a synthetic demo movement track
#'
#' Creates a synthetic `sf` POINT object with mixed migration and
#' stationary behaviour for one or more simulated individuals. Each
#' individual visits **two stationary sites** (mimicking wintering and
#' breeding grounds for a migratory bird), connected by a migration
#' trail. Each stationary site is a tight cluster (sd ~500 m),
#' providing the spatial-occupancy structure that [dbscan_habitats()]
#' delineates at a data-derived scale.
#'
#' Per individual, the breakdown is roughly:
#' \itemize{
#'   \item 35% of points: wintering site (tight cluster at origin).
#'   \item 35% of points: breeding site (tight cluster, +5 deg lon,
#'     +8 deg lat from origin).
#'   \item 30% of points: migration trail connecting the two sites
#'     (high-speed, spatially spread).
#' }
#' Used in examples, website articles, and unit tests.
#'
#' @param n_individuals Integer; number of individuals to simulate
#'   (default 2).
#' @param n_points Integer; total number of points per individual
#'   (default 200).
#' @param seed Integer; random seed for reproducibility (default 1).
#'
#' @return An `sf` POINT object in EPSG:4326 with columns:
#'   \itemize{
#'     \item `individual`: character (e.g., `"ind_01"`).
#'     \item `time`: POSIXct, hourly timestamps starting 2023-01-01.
#'     \item `speed`: numeric, m/s.
#'     \item `lon`, `lat`: numeric coordinates (also retained as columns).
#'   }
#'
#' @examples
#' demo <- make_demo_track()
#' plot(demo["speed"])
#'
#' @export
make_demo_track <- function(n_individuals = 2L,
                            n_points = 200L,
                            seed = 1L) {

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required.", call. = FALSE)
  }

  set.seed(seed)
  out_list <- vector("list", n_individuals)

  for (i in seq_len(n_individuals)) {
    n_migr     <- floor(n_points * 0.30)
    n_stat     <- n_points - n_migr
    n_winter   <- floor(n_stat / 2)
    n_breeding <- n_stat - n_winter

    # Origin offset per individual so multiple individuals are spatially
    # separable when pooled.
    origin_lon <- (i - 1) * 20
    origin_lat <- (i - 1) * 10

    # Wintering site (southern, tight cluster)
    winter_lon <- stats::rnorm(n_winter, mean = origin_lon, sd = 0.005)
    winter_lat <- stats::rnorm(n_winter, mean = origin_lat, sd = 0.005)

    # Breeding site (5 deg E, 8 deg N from origin, tight cluster)
    breeding_lon <- stats::rnorm(n_breeding,
                                 mean = origin_lon + 5, sd = 0.005)
    breeding_lat <- stats::rnorm(n_breeding,
                                 mean = origin_lat + 8, sd = 0.005)

    stat_speed <- stats::rgamma(n_stat, shape = 2, rate = 4)

    # Migration trail: from winter to breeding, fast and spread
    migr_lon <- seq(origin_lon, origin_lon + 5,
                    length.out = n_migr) +
                stats::rnorm(n_migr, sd = 0.02)
    migr_lat <- seq(origin_lat, origin_lat + 8,
                    length.out = n_migr) +
                stats::rnorm(n_migr, sd = 0.02)
    migr_speed <- stats::rgamma(n_migr, shape = 5, rate = 0.3)

    lon   <- c(winter_lon, breeding_lon, migr_lon)
    lat   <- c(winter_lat, breeding_lat, migr_lat)
    speed <- c(stat_speed, migr_speed)
    time  <- seq.POSIXt(as.POSIXct("2023-01-01", tz = "UTC"),
                        by = "hour",
                        length.out = length(lon))

    # Random temporal order (mimic real-world non-sequential tracks)
    ord <- sample.int(length(lon))

    df <- data.frame(
      individual = sprintf("ind_%02d", i),
      time       = time[ord],
      speed      = speed[ord],
      lon        = lon[ord],
      lat        = lat[ord],
      stringsAsFactors = FALSE
    )
    out_list[[i]] <- df
  }

  combined <- do.call(rbind, out_list)
  sf::st_as_sf(combined,
               coords = c("lon", "lat"),
               crs    = 4326,
               remove = FALSE)
}

