test_that("detect_nests finds a tight nest cluster", {
  # UTM-style coords in METERS
  set.seed(1)
  nest_x   <- rnorm(80, mean = 0, sd = 5)     # 5 m sd around nest centre
  nest_y   <- rnorm(80, mean = 0, sd = 5)
  wander_x <- rnorm(20, mean = 0, sd = 500)   # 500 m wandering scale
  wander_y <- rnorm(20, mean = 0, sd = 500)

  pts <- sf::st_as_sf(
    data.frame(
      lon  = c(nest_x, wander_x),
      lat  = c(nest_y, wander_y),
      time = seq.POSIXt(as.POSIXct("2023-06-01"), by = "2 hour",
                        length.out = 100),
      individual = "bird_A"
    ),
    coords = c("lon", "lat"), crs = 32633   # UTM33N (metric)
  )

  res <- detect_nests(pts, time_col = "time",
                     individual_col = "individual",
                     min_pts = 30L, search_dist = 50,
                     time_window = "3 weeks",
                     verbose = FALSE)

  expect_type(res, "list")
  expect_true("points" %in% names(res))
  expect_true("nests"  %in% names(res))
  expect_gte(nrow(res$nests), 1L)
  expect_true(all(c("nest_id", "individual", "start_time",
                    "end_time") %in% names(res$nests)))
})


test_that("detect_nests respects the three-tier intersection", {
  set.seed(2)
  # Dense core (5 m sd) + diffuse halo (200 m sd) in UTM metres
  core_x <- rnorm(100, sd = 5)
  core_y <- rnorm(100, sd = 5)
  halo_x <- rnorm(50,  sd = 200)
  halo_y <- rnorm(50,  sd = 200)

  pts <- sf::st_as_sf(
    data.frame(
      lon = c(core_x, halo_x),
      lat = c(core_y, halo_y),
      time = seq.POSIXt(as.POSIXct("2023-06-01"), by = "1 hour",
                        length.out = 150),
      individual = "bird_B"
    ),
    coords = c("lon", "lat"), crs = 32633
  )

  res <- detect_nests(pts, time_col = "time",
                     individual_col = "individual",
                     min_pts = 40L, search_dist = 60,
                     time_window = "3 weeks",
                     verbose = FALSE)

  # Core points (100, within ~15 m) should survive all 3 tiers.
  # Most halo points (200 m sd) are too far to meet the strictest tier.
  n_nest_pts <- sum(!is.na(res$points$nest_id))
  expect_gte(n_nest_pts, 50L)
  expect_lte(n_nest_pts, 130L)
})


test_that("detect_nests reports no nest when density too low", {
  set.seed(3)
  # 60 points diffusely scattered (500 m sd) -- no nest signal
  pts <- sf::st_as_sf(
    data.frame(
      lon = rnorm(60, sd = 500),
      lat = rnorm(60, sd = 500),
      time = seq.POSIXt(as.POSIXct("2023-06-01"), by = "1 hour",
                        length.out = 60),
      individual = "bird_C"
    ),
    coords = c("lon", "lat"), crs = 32633
  )

  res <- detect_nests(pts, time_col = "time",
                     individual_col = "individual",
                     min_pts = 40L, search_dist = 50,
                     verbose = FALSE)

  expect_equal(nrow(res$nests), 0L)
  expect_true(all(is.na(res$points$nest_id)))
})


test_that("detect_nests rejects invalid inputs", {
  expect_error(detect_nests(data.frame(x = 1)),
               "must be an `sf` object")

  pts <- sf::st_as_sf(data.frame(x = 1, y = 1, t = 1),
                     coords = c("x", "y"), crs = 32633)
  expect_error(detect_nests(pts, time_col = "t"),
               "POSIXct or Date")
})


test_that(".parse_time_window handles common units", {
  expect_equal(.parse_time_window("1 day"),     86400)
  expect_equal(.parse_time_window("3 weeks"),   86400 * 7 * 3)
  expect_equal(.parse_time_window("1 hour"),    3600)
  expect_equal(.parse_time_window("30 minutes"), 1800)
  # Bad unit -> "Unknown time unit"
  expect_error(.parse_time_window("3 fortnights"), "Unknown time unit")
  # Bad number -> "Cannot parse" (use a 2-token string so length check passes)
  expect_error(.parse_time_window("abc weeks"), "Cannot parse")
  # Wrong format -> "must look like"
  expect_error(.parse_time_window("not a window"), "must look like")
})


test_that(".st_dbscan_core matches dbscan::dbscan when time is uniform", {
  set.seed(4)
  coords <- cbind(
    x = rnorm(40, sd = 0.5),
    y = rnorm(40, sd = 0.5)
  )
  times <- rep(as.POSIXct("2023-06-01"), 40)

  res_st <- .st_dbscan_core(coords, times,
                            eps_spatial = 1,
                            eps_temporal = 1e9,
                            min_pts = 5L)
  res_db <- dbscan::dbscan(coords, eps = 1, minPts = 5)

  expect_equal(res_st$n_clusters,
               length(unique(res_db$cluster)) -
                 as.integer(any(res_db$cluster == 0L)))
})

