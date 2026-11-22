test_that("weighted relaxation matches weighted least squares", {
  x <- cbind(seq_len(10), rep(c(1, 0), 5))
  y <- 3 + 2 * x[, 1] - x[, 2] + sin(x[, 1])
  y[1] <- 1000
  weights <- c(0, 1, 2, 3, 2, 1, 4, 5, 1, 2)
  expected <- stats::lm.wfit(cbind(1, x), y, weights)$coefficients

  for (design in list(x, Matrix::Matrix(x, sparse = TRUE))) {
    fit <- SLOPE(
      design,
      y,
      alpha = 0.01,
      gamma = 0,
      weights = weights,
      tol = 1e-9
    )
    rescaled <- SLOPE(
      design,
      y,
      alpha = 0.01,
      gamma = 0,
      weights = 10 * weights,
      tol = 1e-9
    )

    expect_equal(as.numeric(coef(fit)), unname(expected), tolerance = 1e-4)
    expect_equal(coef(fit), coef(rescaled), tolerance = 1e-4)
    expect_identical(stats::weights(fit), weights)
  }
})

test_that("weighted AUC uses the weights of both observations in each pair", {
  x <- matrix(1:8)
  y <- factor(c("no", "no", "yes", "no", "yes", "yes", "no", "yes"))
  weights <- c(0, 1, 4, 2, 1, 5, 1, 3)
  fit <- SLOPE(x, y, family = "binomial", alpha = 0.05, weights = weights)

  prob <- predict(fit, x, type = "response")
  positive <- y == "yes"
  differences <- outer(prob[positive], prob[!positive], `-`)
  pair_weights <- outer(weights[positive], weights[!positive])
  expected <- sum(
    pair_weights * ((differences > 0) + 0.5 * (differences == 0))
  ) /
    sum(pair_weights)

  expect_equal(score(fit, x, y, measure = "auc", weights = weights), expected)
})

test_that("observation weights are validated", {
  x <- matrix(1:8, 4)
  y <- c(1, 2, 4, 3)

  expect_error(SLOPE(x, y, weights = numeric(0)), "one entry per row")
  expect_error(SLOPE(x, y, weights = c(1, 1, 1, NA)), "finite and nonnegative")
  expect_error(SLOPE(x, y, weights = c(1, 1, 1, Inf)), "finite and nonnegative")
  expect_error(SLOPE(x, y, weights = c(1, 1, 1, -1)), "finite and nonnegative")
  expect_error(SLOPE(x, y, weights = rep(0, 4)), "positive total weight")
  expect_error(
    cvSLOPE(x, y, weights = c(1, 0, 0, 0), n_folds = 2),
    "fold [0-9]+ in repeat 1 has zero total observation weight"
  )
  expect_error(
    SLOPE(x, y, weights = 1:4, alpha = "estimate"),
    "Nonuniform observation weights"
  )
})

test_that("binomial fits with an intercept need positive weight in both classes", {
  x <- matrix(1:8)
  y <- rep(c(0, 1), 4)
  expect_error(
    SLOPE(x, y, family = "binomial", weights = 1 - y, alpha = 0.1),
    "positive total weight in both classes"
  )
})

test_that("cross-validation weights the fits and held-out errors", {
  x <- cbind(1:18, cos(1:18))
  y <- 3 + 2 * x[, 1] - x[, 2] + sin(x[, 1])
  weights <- rep(c(0, 1, 3), 6)

  set.seed(42)
  folds <- createFolds(nrow(x), 3)[[1]]
  set.seed(42)
  tune <- cvSLOPE(
    x,
    y,
    n_folds = 3,
    path_length = 3,
    gamma = 1,
    alpha_scale = "relative",
    weights = weights,
    tol = 1e-9
  )
  expected <- lapply(folds, function(fold) {
    test <- fold + 1L
    fit <- SLOPE(
      x[-test, ],
      y[-test],
      q = 0.2,
      alpha = tune$data[[1]]$alphas,
      weights = weights[-test],
      tol = 1e-9
    )
    errors <- (predict(fit, x[test, ]) - y[test])^2
    apply(errors, 2, stats::weighted.mean, w = weights[test])
  })

  expect_equal(tune$data[[1]]$score, do.call(rbind, expected), tolerance = 1e-6)
  expect_equal(tune$model$alpha, tune$optima$alpha)
  expect_identical(stats::weights(tune$model), weights)
  expect_equal(
    coef(refit(tune, tol = 1e-9)),
    coef(tune$model),
    tolerance = 1e-6
  )
  expect_null(stats::weights(refit(tune, weights = NULL)))
})

test_that("legacy training reuses weighted absolute penalties in each fold", {
  x <- cbind(1:18, cos(1:18))
  y <- 3 + 2 * x[, 1] - x[, 2] + sin(x[, 1])
  weights <- rep(c(0, 1, 3), 6)

  set.seed(42)
  tune <- trainSLOPE(
    x,
    y,
    number = 3,
    measure = "mse",
    lambda = "lasso",
    alpha = c(0.6, 0.2),
    alpha_scale = "relative",
    weights = weights,
    tol = 1e-9
  )
  set.seed(42)
  absolute <- trainSLOPE(
    x,
    y,
    number = 3,
    measure = "mse",
    lambda = "lasso",
    alpha = tune$model$alpha,
    weights = weights,
    tol = 1e-9
  )

  expect_equal(tune$summary, absolute$summary, tolerance = 1e-6)
  expect_identical(stats::weights(refit(tune)), weights)
})
