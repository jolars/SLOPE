test_that("regularization path correctly stops if max_variables reached", {

  x <- scale(heart$x)
  y <- heart$y

  fit <- SLOPE(x, y,
               family = "binomial",
               max_variables = 10,
               intercept = FALSE,
               lambda = "bh")

  n_sigma <- length(fit$sigma)

  n_var <- sum(unique(abs(signif(coef(fit)[, n_sigma - 1])), 4) != 0)

  expect_lte(n_var, 10)
})

