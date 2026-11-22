test_that("cross-validation for gaussian and binomial", {
  set.seed(48)

  for (family in c("gaussian", "binomial")) {
    xy <- SLOPE:::randomProblem(1e3, 2, response = family)

    fit <- cvSLOPE(
      xy$x,
      xy$y,
      family = family,
      measure = "mse",
      n_folds = 2,
      q = c(0.1, 0.2),
      n_repeats = 2,
      path_length = 2
    )
    expect_s3_class(fit, "TrainedSLOPE")
  }
})

test_that("cvSLOPE works properly for binomial family", {
  set.seed(1454)
  xy <- SLOPE:::randomProblem(200, p = 10, q = 0.5, response = "binomial")
  x <- xy$x
  y <- xy$y

  fit <- cvSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    n_folds = 2,
    family = "binomial",
    path_length = 20
  )

  expect_true(!anyNA(fit$data$mean))
})

test_that("cvSLOPE works properly for gaussian family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 10, q = 0.5, response = "gaussian")
  x <- xy$x
  y <- xy$y

  fit <- cvSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    n_folds = 2,
    family = "gaussian",
    path_length = 20
  )

  expect_true(!anyNA(fit$optima))
})

test_that("cvSLOPE works properly for poisson family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 5, q = 0.5, response = "poisson")
  x <- xy$x
  y <- xy$y

  fit <- cvSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    n_folds = 2,
    family = "poisson",
    path_length = 20
  )

  expect_true(!anyNA(fit$optima))
})

test_that("cvSLOPE works properly for multinomial family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 5, q = 0.5, response = "multinomial")
  x <- xy$x
  y <- xy$y

  fit <- cvSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    n_folds = 2,
    measure = "auc",
    family = "multinomial",
    path_length = 20
  )

  expect_true(!anyNA(fit$optima))
})

test_that("cvSLOPE returns error in the case of invalid measures", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 100, q = 0.5, response = "gaussian")
  x <- xy$x
  y <- xy$y

  expect_error(cvSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    measure = "misclass",
    family = "gaussian",
    path_length = 20
  ))
})

# TODO: Find out why this test fails
# test_that("Returned AUC from cvSLOPE is correct", {
#   set.seed(42)
#   xy <- SLOPE:::randomProblem(200, p = 100, q = 0.5, response = "binomial")
#   x <- xy$x
#   y <- xy$y
#
#   tuned <- cvSLOPE(
#     x,
#     y,
#     q = c(0.1, 0.2),
#     measure = "auc",
#     family = "binomial",
#   )
#
#   expect_gte(tuned$optima$mean[1], 0)
#   expect_lte(tuned$optima$mean[1], 1)
# })

test_that("Cross-validating on gamma works", {
  set.seed(35)

  xy <- SLOPE:::randomProblem(200, p = 2)
  x <- xy$x
  y <- xy$y

  tuned <- cvSLOPE(x, y, q = c(0.1, 0.2), gamma = c(0.0, 1.0), )

  expect_gte(tuned$optima$mean[1], 0)
  expect_lte(tuned$optima$mean[1], 1)
})

test_that("cross-validation uses supplied alpha candidates and refits the winner", {
  x <- cbind(1:18, cos(1:18))
  y <- 3 + 2 * x[, 1] - x[, 2] + sin(x[, 1])
  weights <- rep(c(0, 1, 3), 6)
  candidates <- c(0.6, 0.2, 0)

  for (design in list(x, Matrix::Matrix(x, sparse = TRUE))) {
    set.seed(42)
    folds <- createFolds(nrow(x), 3)[[1]]
    set.seed(42)
    tune <- cvSLOPE(
      design,
      y,
      q = c(0.1, 0.2),
      alpha = candidates,
      alpha_scale = "relative",
      n_folds = 3,
      gamma = 1,
      weights = weights,
      tol = 1e-9
    )

    for (grid in tune$data) {
      q <- grid$params[["q"]]
      full_fit <- SLOPE(
        design,
        y,
        q = q,
        alpha = candidates,
        alpha_scale = "relative",
        weights = weights,
        tol = 1e-9
      )
      expect_equal(grid$alphas, full_fit$alpha)

      expected <- lapply(folds, function(fold) {
        test <- fold + 1L
        fit <- SLOPE(
          design[-test, , drop = FALSE],
          y[-test],
          q = q,
          alpha = full_fit$alpha,
          weights = weights[-test],
          tol = 1e-9
        )
        errors <- (predict(fit, design[test, , drop = FALSE]) - y[test])^2
        apply(errors, 2, stats::weighted.mean, w = weights[test])
      })
      expect_equal(grid$score, do.call(rbind, expected), tolerance = 1e-6)
    }

    expect_equal(tune$model$alpha, tune$optima$alpha)
  }

  absolute <- cvSLOPE(x, y, alpha = 0.123, n_folds = 3, gamma = 1)
  expect_equal(absolute$data[[1]]$alphas, 0.123)
  expect_equal(absolute$model$alpha, 0.123)
})

test_that("relative CV alpha uses the supplied lambda sequence", {
  x <- cbind(1:18, cos(1:18))
  y <- 3 + 2 * x[, 1] - x[, 2] + sin(x[, 1])
  args <- list(x = x, y = y, alpha = c(0.6, 0.2), lambda = c(4, 2), gamma = 1)
  fit <- do.call(SLOPE, c(args, list(alpha_scale = "relative")))

  set.seed(42)
  relative <- do.call(cvSLOPE, c(args, list(alpha_scale = "relative")))
  set.seed(42)
  args$alpha <- fit$alpha
  absolute <- do.call(cvSLOPE, args)

  expect_equal(relative$data[[1]]$alphas, fit$alpha)
  expect_equal(relative$summary, absolute$summary)
  expect_equal(coef(relative$model), coef(absolute$model))
})
