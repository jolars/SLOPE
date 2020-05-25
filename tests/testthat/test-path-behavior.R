test_that("regularization path correctly stops if max_variables reached", {

  x <- scale(heart$x)
  y <- heart$y

  fit <- SLOPE(x, y,
               family = "binomial",
               max_variables = 10,
               intercept = FALSE,
               lambda = "bh")

  path_length <- length(fit$alpha)

  n_var <- sum(unique(abs(signif(coef(fit)[, path_length - 1])), 4) != 0)

  expect_lte(n_var, 10)
})

