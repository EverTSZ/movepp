test_that("detect_dominant_axis returns expected structure", {
  demo <- make_demo_track(n_individuals = 1, n_points = 200, seed = 1L)
  result <- detect_dominant_axis(demo)
  expect_s3_class(result, "movepp_dominant_axis")
  expect_true(result$primary_axis %in% c("latitude", "longitude"))
  expect_true(result$pc1_variance_explained >= 0 &&
              result$pc1_variance_explained <= 1)
})

test_that("detect_dominant_axis identifies latitude for bird demo", {
  demo <- make_demo_track(n_individuals = 1, n_points = 300, seed = 1L)
  result <- detect_dominant_axis(demo)
  expect_equal(result$primary_axis, "latitude")
})

test_that("detect_dominant_axis rejects non-sf input", {
  expect_error(detect_dominant_axis(data.frame(x = 1:5, y = 1:5)),
               "must be an `sf` object")
})

test_that("compute_mpi requires a time column", {
  demo <- make_demo_track(n_individuals = 1, n_points = 200, seed = 1L)
  expect_error(compute_mpi(demo, time_col = NULL, verbose = FALSE))
})

test_that("compute_mpi adds an mpi column in [0,1]", {
  demo <- make_demo_track(n_individuals = 1, n_points = 300, seed = 1L)
  r <- compute_mpi(demo, time_col = "time", verbose = FALSE)
  expect_true("mpi" %in% names(r))
  expect_true(all(r$mpi >= 0 & r$mpi <= 1, na.rm = TRUE))
})

test_that("compute_mpi accepts supplied population anchors", {
  demo <- make_demo_track(n_individuals = 1, n_points = 300, seed = 1L)
  r <- compute_mpi(demo, time_col = "time",
                   winter_doy = 15, breed_doy = 175, verbose = FALSE)
  expect_equal(attr(r, "mpi_winter_doy"), 15)
  expect_equal(attr(r, "mpi_breed_doy"), 175)
  expect_true(all(r$mpi >= 0 & r$mpi <= 1, na.rm = TRUE))
})
