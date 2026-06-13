test_that("balm_segmentation classifies movement states", {
  demo <- make_demo_track(n_individuals = 2, n_points = 400, seed = 1L)
  g    <- compute_step_speed(demo, time_col = "time", individual_col = "individual")
  res  <- balm_segmentation(g, variable_col = "step_speed",
                            individual_col = "individual", verbose = FALSE)
  expect_s3_class(res, "sf")
  expect_true(all(c("balm_deviation", "balm_lag", "cluster_type") %in% names(res)))
  expect_true(is.factor(res$cluster_type))
  expect_true(all(levels(res$cluster_type) %in% c("HH", "LL", "HL", "LH")))
  expect_false(is.null(attr(res, "balm_reference")))
})
