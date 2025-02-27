test_that("diagnostics are working properly", {
  set.seed(33)
  xy <- SLOPE:::randomProblem(100, 2, q = 1)

  fit <- SLOPE(xy$x, xy$y, diagnostics = TRUE, alpha = 1)

  expect_is(fit$diagnostics, "data.frame")
  tmp <- tempfile()
  grDevices::png(tmp)
  expect_silent(plotDiagnostics(fit, xvar = "iteration"))
  grDevices::dev.off()
  unlink(tmp)
})
