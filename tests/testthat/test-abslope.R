
test_that("ABSLOPE() returns proper output.", {

  set.seed(17)
  xy <- SLOPE:::randomProblem(1e2, 200, response = "gaussian")
  X <- as.matrix(xy$x)
  Y <- xy$y
  fit <- ABSLOPE(X, Y)
  expect_equal(sum(fit$beta != 0), 6)
  expect_equal(fit$sigma, 21.6073367639117)
  expect_equal(fit$intercept, -0.0135406689752608)
})
