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

test_that("TrainedSLOPE summary works correctly", {
  tune <- cvSLOPE(
    subset(mtcars, select = c("mpg", "drat", "wt")),
    mtcars$hp,
    q = c(0.1, 0.2),
    n_folds = 5,
    n_repeats = 2
  )

  s <- summary(tune)

  expect_s3_class(s, "summary_TrainedSLOPE")
  expect_true("call" %in% names(s))
  expect_true("measure" %in% names(s))
  expect_true("optima" %in% names(s))
  expect_true("n_folds" %in% names(s))
  expect_true("n_repeats" %in% names(s))
  expect_true("n_models" %in% names(s))

  expect_equal(s$n_folds, 5)
  expect_equal(s$n_repeats, 2)
})

test_that("print.summary_TrainedSLOPE runs without error", {
  tune <- cvSLOPE(
    subset(mtcars, select = c("mpg", "drat", "wt")),
    mtcars$hp,
    q = c(0.1, 0.2),
    n_folds = 5
  )

  s <- summary(tune)

  expect_output(print(s), "Cross-validation:")
  expect_output(print(s), "Folds:")
  expect_output(print(s), "Models evaluated:")
  expect_output(print(s), "Performance measure:")
  expect_output(print(s), "Optimal parameters:")
})

test_that("refit works for TrainedSLOPE objects", {
  x <- bodyfat$x
  y <- bodyfat$y
  
  tune <- trainSLOPE(x, y, q = c(0.1, 0.2), measure = c("mse", "mae"), number = 3)
  
  fit1 <- refit(tune, x, y)
  expect_s3_class(fit1, "SLOPE")
  expect_s3_class(fit1, "GaussianSLOPE")
  
  expect_equal(length(fit1$alpha), 1)
  
  fit2 <- refit(tune, x, y, measure = "mae")
  expect_s3_class(fit2, "SLOPE")
})

test_that("refit errors on invalid measure", {
  tune <- trainSLOPE(bodyfat$x, bodyfat$y, q = 0.1, measure = "mse", number = 3)
  
  expect_error(
    refit(tune, bodyfat$x, bodyfat$y, measure = "invalid"),
    "measure 'invalid' not found"
  )
})

test_that("cvSLOPE stores refitted model when refit = TRUE", {
  tune <- cvSLOPE(
    bodyfat$x,
    bodyfat$y,
    q = c(0.1, 0.2),
    n_folds = 3,
    refit = TRUE
  )
  
  expect_true("model" %in% names(tune))
  expect_s3_class(tune$model, "SLOPE")
})

test_that("cvSLOPE does not store model when refit = FALSE", {
  tune <- cvSLOPE(
    bodyfat$x,
    bodyfat$y,
    q = c(0.1, 0.2),
    n_folds = 3,
    refit = FALSE
  )
  
  expect_false("model" %in% names(tune))
})

test_that("trainSLOPE still stores model", {
  tune <- trainSLOPE(bodyfat$x, bodyfat$y, q = 0.1, number = 3)
  
  expect_true("model" %in% names(tune))
  expect_s3_class(tune$model, "SLOPE")
})
