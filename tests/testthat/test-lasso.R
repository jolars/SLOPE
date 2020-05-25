test_that("lasso and slope fits are equivalent if all lambda are equal", {
  library(glmnet)
  glmnet.control(fdev = 0)

  set.seed(1)
  xy <- SLOPE:::randomProblem(100, 10)
  x <- xy$x
  y <- xy$y

  lasso <- glmnet(x, y, standardize = FALSE)
  lambda <- lasso$lambda

  slope <- SLOPE(x, y,
                 center = FALSE,
                 scale = FALSE,
                 lambda = rep(lambda[1], ncol(x)),
                 alpha = exp(seq(log(1), log(1e-4), length.out = 100)))

  lasso_coef <- coef(lasso)
  slope_coef <- coef(slope)

  expect_equivalent(as.matrix(lasso_coef)[, 1:75],
                    as.matrix(slope_coef)[, 1:75], tol = 1e-3)
})
