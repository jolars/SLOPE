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

  expect_error(
    cvSLOPE(
      x,
      y,
      q = c(0.1, 0.2),
      measure = "misclass",
      family = "gaussian",
      path_length = 20
    )
  )
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

  tuned <- cvSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    gamma = c(0.0, 1.0),
  )

  expect_gte(tuned$optima$mean[1], 0)
  expect_lte(tuned$optima$mean[1], 1)
})
