test_that("lasso and slope fits are equivalent if all lambda are equal", {

  set.seed(1)
  xy <- SLOPE:::randomProblem(100, 10)
  x <- xy$x
  y <- xy$y

  lambda <- 0.6

  library(glmnet)
  glmnet.control(fdev = 0)
  lasso <- glmnet(x, y, lambda = 0.6, standardize = FALSE)
  lasso_coef <- as.vector(coef(lasso))

  slope <- SLOPE(x, y,
                 center = FALSE,
                 scale = FALSE,
                 lambda = rep(lambda, ncol(x)),
                 alpha = 1)

  slope_coef <- coef(slope)

  expect_equivalent(lasso_coef, slope_coef, tol = 1e-3)
})
