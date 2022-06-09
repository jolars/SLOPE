test_that("diagnostics are working properly", {
  set.seed(33)
  xy <- SLOPE:::randomProblem(100, 2, q = 1)

  fit <- SLOPE(xy$x, xy$y, diagnostics = TRUE, path_length = 10, alpha = 1)

  expect_is(fit$diagnostics, "data.frame")
  p <- plotDiagnostics(fit, xvar = "iteration")
  expect_s3_class(p, "ggplot")

  skip_on_ci()

  vdiffr::expect_doppelganger("plotDiagnostics-in-test", p)
})
