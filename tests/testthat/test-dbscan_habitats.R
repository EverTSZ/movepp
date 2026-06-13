test_that("dbscan_habitats requires eps and minPts (no defaults)", {
  demo <- make_demo_track(n_individuals = 2, n_points = 400)
  seg  <- balm_segmentation(demo, "speed", "individual", verbose = FALSE)
  stat <- seg[!is.na(seg$cluster_type) & seg$cluster_type %in% c("LL","LH"), ]
  expect_error(dbscan_habitats(stat, individual_col = "individual"))
})

test_that("dbscan_habitats returns an sf with a per-individual cluster_id", {
  demo <- make_demo_track(n_individuals = 2, n_points = 400)
  seg  <- balm_segmentation(demo, "speed", "individual", verbose = FALSE)
  stat <- seg[!is.na(seg$cluster_type) & seg$cluster_type %in% c("LL","LH"), ]
  hab <- dbscan_habitats(stat, individual_col = "individual", eps = 10, minPts = 3)
  expect_s3_class(hab, "sf")
  expect_true("cluster_id" %in% names(hab))
  expect_true(all(hab$cluster_id >= 1L))
})

test_that("detect_habitat_params needs a track (or eps) to derive eps", {
  demo <- make_demo_track(n_individuals = 1, n_points = 120)
  expect_error(detect_habitat_params(demo, individual_col = "individual"))
})

test_that("detect_habitat_params derives eps and minPts on demo data", {
  g <- compute_step_speed(godwit_demo, time_col = "time",
                          individual_col = "individual", direction = "centered")
  seg <- balm_segmentation(g, "step_speed", "individual", verbose = FALSE)
  stat <- seg[!is.na(seg$cluster_type) & seg$cluster_type %in% c("LL","LH"), ]
  p <- detect_habitat_params(stat, track = seg, individual_col = "individual",
                             time_col = "time", verbose = FALSE)
  expect_true(is.numeric(p$eps) && p$eps > 0)
  expect_true(is.numeric(p$minPts) && p$minPts >= 2)
})
