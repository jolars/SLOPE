context("summary")

test_that("SLOPE object contains n_observations and n_predictors", {
  fit <- SLOPE(bodyfat$x, bodyfat$y)

  expect_true("n_observations" %in% names(fit))
  expect_true("n_predictors" %in% names(fit))
  expect_equal(fit$n_observations, nrow(bodyfat$x))
  expect_equal(fit$n_predictors, ncol(bodyfat$x))
})

test_that("summary displays observations and predictors correctly", {
  fit <- SLOPE(heart$x, heart$y)
  s <- summary(fit)

  expect_equal(s$n_obs, nrow(heart$x))
  expect_equal(s$n_predictors, ncol(heart$x))
})

test_that("summary works for different families", {
  # Gaussian
  fit_gaussian <- SLOPE(bodyfat$x, bodyfat$y)
  s_gaussian <- summary(fit_gaussian)
  expect_equal(s_gaussian$n_obs, nrow(bodyfat$x))
  expect_equal(s_gaussian$n_predictors, ncol(bodyfat$x))

  # Binomial
  fit_binomial <- SLOPE(heart$x, heart$y, family = "binomial")
  s_binomial <- summary(fit_binomial)
  expect_equal(s_binomial$n_obs, nrow(heart$x))
  expect_equal(s_binomial$n_predictors, ncol(heart$x))

  # Multinomial
  fit_multi <- SLOPE(wine$x, wine$y, family = "multinomial")
  s_multi <- summary(fit_multi)
  expect_equal(s_multi$n_obs, nrow(wine$x))
  expect_equal(s_multi$n_predictors, ncol(wine$x))
})

test_that("print.summary_SLOPE runs without error", {
  fit <- SLOPE(heart$x, heart$y)
  s <- summary(fit)

  expect_output(print(s), "Observations:")
  expect_output(print(s), "Predictors:")
  expect_output(print(s), "Path summary")
})
