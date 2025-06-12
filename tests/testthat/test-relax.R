test_that("relaxed SLOPE", {
  set.seed(30)
  xy <- SLOPE:::randomProblem(100, 10)

  x <- diag(3)
  beta <- c(1, 0.0, 1)
  y <- x %*% beta

  # one parameter
  fit <- SLOPE(xy$x, xy$y, gamma = 0.0001, scale = FALSE, center = FALSE)
  expect_silent(dont_plot(fit))
})
