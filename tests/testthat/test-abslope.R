
test_that("ABSLOPE() returns proper output.", {

  set.seed(17)
  xy <- SLOPE:::randomProblem(1e2, 200, response = "gaussian")
  X <- as.matrix(xy$x)
  Y <- xy$y
  fit <- ABSLOPE(X, Y)
  expect_equal(which(fit$coefficients != 0), c(`(Intercept)` = 1L, V57 = 58L,
                                               V61 = 62L, V96 = 97L, V99 = 100L,
                                               V126 = 127L, V188 = 189L))
  expect_equal(fit$sigma, 21.64159, tolerance = 4)
})


