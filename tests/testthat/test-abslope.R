

test_that("ABSLOPE", {

  set.seed(17)
  xy <- SLOPE:::randomProblem(1e2, 2, response = "gaussian")
  X <- as.matrix(xy$x)
  Y <- xy$y
  fit <- ABSLOPE(X, Y)

  expect_false(all(is.na(fit$beta)))
})
