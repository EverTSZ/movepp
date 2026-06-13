# tests/testthat/test-classify_phases.R
# Replaces the old test-gmm_voting.R (classify_phases_gmm / test_phase_normality
# were removed). Self-contained: builds a synthetic "spike-plateau-spike" MPI
# rather than depending on the full upstream pipeline.

test_that("classify_phases cuts a 3-mode MPI into wintering/stopover/breeding", {
  skip_if_not_installed("sf")
  set.seed(1)
  # wintering ~0, broad stopover plateau, tight breeding spike ~1
  mpi <- c(stats::rnorm(120, 0.03, 0.02),
           stats::runif( 80, 0.30, 0.70),
           stats::rnorm( 60, 0.97, 0.02))
  mpi <- pmin(1, pmax(0, mpi))
  n   <- length(mpi)
  pts <- sf::st_as_sf(
    data.frame(mpi        = mpi,
               individual = "A",
               cluster_id = rep(seq_len(ceiling(n / 5)), each = 5)[seq_len(n)],
               x = stats::runif(n), y = stats::runif(n)),
    coords = c("x", "y"), crs = 4326)

  res <- classify_phases(pts, n_phases = 3, verbose = FALSE)

  expect_s3_class(res, "sf")
  expect_true(all(c("phase_label", "phase_voted_label", "phase") %in% names(res)))

  fit <- attr(res, "phase_fit")
  expect_equal(fit$n_phases, 3L)
  expect_length(fit$phase_cuts, 2L)
  expect_setequal(levels(res$phase), c("Wintering", "Stopover", "Breeding"))

  # points near 1 must land in Breeding; points near 0 in Wintering
  expect_gt(mean(res$phase[mpi > 0.9] == "Breeding"),  0.8)
  expect_gt(mean(res$phase[mpi < 0.1] == "Wintering"), 0.8)
})

test_that("classify_phases degrades to 2 phases when only one valley exists", {
  skip_if_not_installed("sf")
  set.seed(2)
  mpi <- c(stats::rnorm(100, 0.05, 0.03), stats::rnorm(100, 0.95, 0.03))
  mpi <- pmin(1, pmax(0, mpi))
  n   <- length(mpi)
  pts <- sf::st_as_sf(
    data.frame(mpi = mpi, individual = "A",
               cluster_id = rep(seq_len(40), length.out = n),
               x = stats::runif(n), y = stats::runif(n)),
    coords = c("x", "y"), crs = 4326)

  res <- classify_phases(pts, n_phases = 3, verbose = FALSE)
  expect_lte(attr(res, "phase_fit")$n_phases, 2L)
})

test_that("classify_phases errors on too few MPI values", {
  skip_if_not_installed("sf")
  pts <- sf::st_as_sf(
    data.frame(mpi = runif(5), individual = "A", cluster_id = 1:5,
               x = runif(5), y = runif(5)),
    coords = c("x", "y"), crs = 4326)
  expect_error(classify_phases(pts, verbose = FALSE), "at least 10")
})
