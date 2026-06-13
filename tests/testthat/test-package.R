test_that("package version is reachable", {
  v <- movepp_version()
  expect_s3_class(v, "package_version")
  expect_true(as.character(v) >= "0.1.0")
})
