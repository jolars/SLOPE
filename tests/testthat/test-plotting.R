test_that("plotting works", {
  set.seed(1)
  xy <- SLOPE:::randomProblem(100, 2)

  # one parameter
  fit <- SLOPE(xy$x, xy$y, alpha = 0.2)
  expect_silent(dont_plot(fit))

  # more parameters
  fit <- SLOPE(xy$x, xy$y, path_length = 10)
  expect_silent(dont_plot(fit))
})


test_that("plot.SLOPE works as expected", {
  fit <- SLOPE(heart$x, heart$y)
  p <- plot(fit)
  vdiffr::expect_doppelganger("plot.SLOPE-in-test", p)

  p <- plot(fit, intercept = TRUE, x_variable = "deviance_ratio")
  vdiffr::expect_doppelganger("plot.SLOPE-parameters-in-test", p)
})
