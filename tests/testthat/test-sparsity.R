test_that("sparse and dense implementations give equivalent results", {
  set.seed(2)
  n <- 100
  p <- 2

  for (family in c("gaussian", "binomial", "poisson")) {
    d <- SLOPE:::randomProblem(n, p, 0.5, density = 0.5, response = family)
    sparse_x <- d$x
    dense_x <- as.matrix(sparse_x)
    y <- d$y
    beta <- d$beta

    sparse_fit <- SLOPE(sparse_x, y, family = family, tol = 1e-8)
    dense_fit <- SLOPE(dense_x, y, family = family, tol = 1e-8)

    sparse_coefs <- as.matrix(coef(sparse_fit))
    dense_coefs <- as.matrix(coef(dense_fit))

    expect_equal(sparse_coefs, dense_coefs, tol = 1e-4)
  }
})
