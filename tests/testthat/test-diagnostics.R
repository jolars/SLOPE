test_that("diagnostics are working properly", {
  xy <- SLOPE:::randomProblem(100, 2, q = 1)

  fit <- SLOPE(xy$x, xy$y, diagnostics = TRUE, n_sigma = 1, sigma = 1)

  expect_is(fit$diagnostics, "data.frame")
  p <- plotDiagnostics(fit)
  expect_s3_class(p, "trellis")
})
