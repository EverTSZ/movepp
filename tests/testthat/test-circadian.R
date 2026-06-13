test_that("compute_solar_elevation adds solar_elevation column", {
  pts <- sf::st_as_sf(data.frame(
    x = c(0, 0), y = c(0, 0),
    time = as.POSIXct(c("2023-06-21 12:00:00", "2023-06-21 00:00:00"),
                       tz = "UTC")
  ), coords = c("x", "y"), crs = 4326)

  res <- compute_solar_elevation(pts, time_col = "time",
                                 verbose = FALSE)
  expect_true("solar_elevation" %in% names(res))
  expect_type(res$solar_elevation, "double")
  # Noon at equator on summer solstice → sun nearly directly overhead
  expect_gt(res$solar_elevation[1], 60)
  # Midnight at equator → sun well below horizon
  expect_lt(res$solar_elevation[2], -60)
})


test_that("classify_circadian assigns expected categories", {
  pts <- sf::st_as_sf(data.frame(
    x = 0, y = 0,
    solar_elevation = c(45, 30, 0, -3, -6, -20)
  ), coords = c("x", "y"), crs = 4326)

  res <- classify_circadian(pts, verbose = FALSE)

  expect_s3_class(res$circadian, "factor")
  expect_equal(levels(res$circadian), c("Day", "Twilight", "Night"))
  expect_equal(as.character(res$circadian),
               c("Day", "Day", "Twilight", "Twilight", "Night", "Night"))
})


test_that("classify_circadian respects custom thresholds", {
  pts <- sf::st_as_sf(data.frame(
    x = 0, y = 0,
    solar_elevation = c(10, -5, -15)
  ), coords = c("x", "y"), crs = 4326)

  res <- classify_circadian(pts, day_threshold = 5,
                            twilight_threshold = -12,
                            verbose = FALSE)
  expect_equal(as.character(res$circadian),
               c("Day", "Twilight", "Night"))
})


test_that("classify_circadian rejects invalid thresholds", {
  pts <- sf::st_as_sf(data.frame(
    x = 0, y = 0, solar_elevation = 0
  ), coords = c("x", "y"), crs = 4326)
  expect_error(classify_circadian(pts, day_threshold = -6,
                                  twilight_threshold = 0),
               "twilight_threshold` must be <")
})


test_that("compute_solar_elevation rejects non-POSIXct time", {
  pts <- sf::st_as_sf(data.frame(
    x = 0, y = 0, time = "2023-06-21"
  ), coords = c("x", "y"), crs = 4326)
  expect_error(compute_solar_elevation(pts, time_col = "time"),
               "POSIXct")
})


test_that("classify_habitat_circadian types habitats by majority period", {
  circ <- c(rep("Night", 8), rep("Day", 2),                 # cluster 1 -> Night
            rep("Day", 5),   rep("Night", 5),               # cluster 2 -> Mixed
            rep("Day", 6),   rep("Twilight", 2), rep("Night", 2))  # 3 -> Day
  cid  <- c(rep(1L, 10), rep(2L, 10), rep(3L, 10))
  n    <- length(circ)
  pts  <- sf::st_as_sf(
    data.frame(individual = "A", cluster_id = cid, circadian = circ,
               x = stats::runif(n), y = stats::runif(n)),
    coords = c("x", "y"), crs = 4326)

  res <- classify_habitat_circadian(pts, individual_col = "individual",
                                    verbose = FALSE)
  expect_true("habitat_circadian" %in% names(res))
  expect_s3_class(res$habitat_circadian, "factor")

  # the type is constant within each habitat
  by_cluster <- tapply(as.character(res$habitat_circadian), res$cluster_id,
                       function(z) unique(z))
  expect_equal(unname(by_cluster[["1"]]), "Night")
  expect_equal(unname(by_cluster[["2"]]), "Mixed")
  expect_equal(unname(by_cluster[["3"]]), "Day")

  summ <- attr(res, "habitat_circadian")
  expect_equal(nrow(summ), 3L)
  expect_true(all(c("cluster_id", "n", "type") %in% names(summ)))
})


test_that("classify_habitat_circadian validates inputs", {
  pts <- sf::st_as_sf(
    data.frame(cluster_id = 1, circadian = "Day", x = 0, y = 0),
    coords = c("x", "y"), crs = 4326)
  expect_error(classify_habitat_circadian(pts, majority = 1.5), "majority")
  expect_error(classify_habitat_circadian(data.frame(a = 1)),
               "must be an `sf`")
})

