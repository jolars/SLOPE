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
