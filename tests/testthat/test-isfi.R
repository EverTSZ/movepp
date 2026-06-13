test_that("compute_isfi computes correct ISFI for known data", {
  # 3 clusters, 3 years
  # Cluster 1 visited in all 3 years
  # Cluster 2 visited in 2 years
  # Cluster 3 visited in 1 year
  # sum_cluster_years = 3 + 2 + 1 = 6
  # n = 3, c = 3
  # ISFI = (6 - 3) / ((3-1) * 3) = 3/6 = 0.5
  df <- data.frame(
    individual = "A",
    cluster_id = c(1, 1, 1, 2, 2, 3),
    time = as.POSIXct(c("2020-01-01", "2021-01-01", "2022-01-01",
                         "2020-06-01", "2021-06-01", "2022-06-01"))
  )
  # st_as_sf with raw lists is fiddly; use explicit columns
  df$x <- 1:6; df$y <- 1:6
  pts <- sf::st_as_sf(df, coords = c("x","y"), crs = 4326)

  res <- compute_isfi(pts, time_col = "time", verbose = FALSE)

  expect_equal(res$n_clusters, 3L)
  expect_equal(res$n_years, 3L)
  expect_equal(res$sum_cluster_years, 6L)
  expect_equal(res$isfi, 0.5, tolerance = 1e-8)
})


test_that("compute_isfi gives 1.0 when every cluster visited every year", {
  df <- data.frame(
    individual = "A",
    cluster_id = c(1, 1, 1, 2, 2, 2),
    time = as.POSIXct(c("2020-01-01", "2021-01-01", "2022-01-01",
                         "2020-06-01", "2021-06-01", "2022-06-01")),
    x = 1:6, y = 1:6
  )
  pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = 4326)
  res <- compute_isfi(pts, time_col = "time", verbose = FALSE)
  expect_equal(res$isfi, 1.0)
})


test_that("compute_isfi gives 0.0 when each cluster only one year", {
  df <- data.frame(
    individual = "A",
    cluster_id = c(1, 2, 3, 4),
    time = as.POSIXct(c("2020-01-01", "2020-02-01",
                         "2021-01-01", "2021-02-01")),
    x = 1:4, y = 1:4
  )
  pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = 4326)
  res <- compute_isfi(pts, time_col = "time", verbose = FALSE)
  # n = 4, c = 2, sum_cluster_years = 4
  # ISFI = (4 - 4) / (1 * 4) = 0
  expect_equal(res$isfi, 0.0)
})


test_that("compute_isfi returns NA when only one year available", {
  df <- data.frame(
    individual = "A",
    cluster_id = c(1, 1, 2),
    time = as.POSIXct(c("2020-01-01", "2020-02-01", "2020-03-01")),
    x = 1:3, y = 1:3
  )
  pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = 4326)
  res <- suppressMessages(
    compute_isfi(pts, time_col = "time", verbose = TRUE)
  )
  expect_true(is.na(res$isfi))
})


test_that("compute_isfi computes per individual when multiple given", {
  df <- data.frame(
    individual = c(rep("A", 4), rep("B", 4)),
    cluster_id = c(1, 1, 2, 2, 1, 2, 3, 4),
    time = rep(as.POSIXct(c("2020-01-01", "2021-01-01")), 4),
    x = 1:8, y = 1:8
  )
  pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = 4326)
  res <- compute_isfi(pts, time_col = "time", verbose = FALSE)
  expect_equal(nrow(res), 2L)
  expect_setequal(res$individual, c("A", "B"))
})


test_that("compute_isfi stratifies by phase when phase_col supplied", {
  df <- data.frame(
    individual = rep("A", 6),
    cluster_id = c(1, 1, 2, 3, 3, 4),
    phase      = rep(c("breeding", "wintering"), each = 3),
    time = as.POSIXct(c("2020-06-01", "2021-06-01", "2022-06-01",
                         "2020-12-01", "2021-12-01", "2022-12-01")),
    x = 1:6, y = 1:6
  )
  pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = 4326)
  res <- compute_isfi(pts, time_col = "time", phase_col = "phase",
                     verbose = FALSE)
  expect_equal(nrow(res), 2L)
  expect_true("phase" %in% names(res))
  expect_setequal(res$phase, c("breeding", "wintering"))
})


test_that("compute_isfi rejects invalid inputs", {
  expect_error(compute_isfi(list(x = 1), time_col = "time"),
               "must be an `sf` or `data.frame`")

  pts <- sf::st_as_sf(
    data.frame(individual = "A", cluster_id = 1, x = 1, y = 1),
    coords = c("x", "y"), crs = 4326)
  expect_error(compute_isfi(pts, time_col = "nope"), "not found")

  pts$bad_time <- "not a time"
  expect_error(compute_isfi(pts, time_col = "bad_time"),
               "POSIXct or Date")
})


test_that("compute_isfi works on data.frame (no sf) for INFI use case", {
  df <- data.frame(
    individual = "A",
    nest_id    = c(1, 1, 2),
    time       = as.POSIXct(c("2020-06-01", "2021-06-01", "2020-07-01"))
  )
  res <- compute_isfi(df, time_col = "time",
                     cluster_col = "nest_id", verbose = FALSE)
  expect_equal(res$n_clusters, 2L)
  expect_equal(res$n_years, 2L)
  expect_equal(res$sum_cluster_years, 3L)
  # ISFI = (3 - 2) / ((2-1) * 2) = 0.5
  expect_equal(res$isfi, 0.5)
})

