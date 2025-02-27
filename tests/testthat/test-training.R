test_that("model training works with trainSLOPE", {
  set.seed(48)

  for (family in c("gaussian", "binomial")) {
    xy <- SLOPE:::randomProblem(1e2, 2, response = family)

    fit <- trainSLOPE(xy$x,
      xy$y,
      family = family,
      number = 2,
      q = c(0.1, 0.2),
      repeats = 2,
      path_length = 2
    )
    expect_s3_class(fit, "TrainedSLOPE")
  }
})


test_that("trainSLOPE works properly for binomial family", {
  set.seed(1454)
  xy <- SLOPE:::randomProblem(200, p = 10, q = 0.5, response = "binomial")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    number = 2,
    family = "binomial",
    path_length = 20
  )

  expect_equal(
    fit$measure$measure,
    c("mse", "mae", "deviance", "misclass", "auc")
  )
  expect_equal(
    fit$measure$label,
    c(
      "Mean Squared Error", "Mean Absolute Error",
      "Binomial Deviance", "Misclassification Rate", "AUC"
    )
  )

  expect_equal(
    fit$optima$mean,
    c(
      0.9811094, 0.361433768605148, 0.179544219401791, 0.07,
      0.112826460692166
    ),
    tolerance = 0.0001
  )
})


test_that("trainSLOPE works properly for gaussian family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 10, q = 0.5, response = "gaussian")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    number = 2,
    family = "gaussian",
    path_length = 20
  )

  expect_equal(
    fit$measure$measure,
    c("mse", "mae")
  )
  expect_equal(
    fit$measure$label,
    c("Mean Squared Error", "Mean Absolute Error")
  )

  expect_equal(
    fit$optima$mean,
    c(0.815915369960497, 1.04529288160428),
    tolerance = 0.001
  )
})

test_that("trainSLOPE works properly for poisson family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 5, q = 0.5, response = "poisson")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    number = 2,
    family = "poisson",
    path_length = 20
  )

  expect_equal(fit$measure$measure, c("mse", "mae"))
  expect_equal(
    fit$measure$label,
    c("Mean Squared Error", "Mean Absolute Error")
  )

  expect_equal(
    fit$optima$mean,
    c(2.17204837325001, 31.7693800794981),
    tolerance = 0.001
  )
})


test_that("trainSLOPE works properly for multinomial family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 5, q = 0.5, response = "multinomial")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    number = 2,
    family = "multinomial",
    path_length = 20
  )

  expect_equal(fit$measure$measure, c("mse", "mae", "deviance", "misclass"))
  expect_equal(fit$measure$label, c(
    "Mean Squared Error",
    "Mean Absolute Error",
    "Multinomial Deviance",
    "Misclassification Rate"
  ))

  expect_equal(
    fit$optima$mean,
    c(55.1784663057681, 0.103107966843569, 0.125, 0.0548346107393888),
    tolerance = 0.0001
  )
})



test_that("trainSLOPE returns error in the case of invalid measures", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 100, q = 0.5, response = "gaussian")
  x <- xy$x
  y <- xy$y

  expect_error(
    trainSLOPE(x, y,
      q = c(0.1, 0.2),
      measure = "misclass",
      family = "gaussian",
      path_length = 20
    ),
    "For the given family: gaussian, measure needs to be one of: mse, mae"
  )
})

test_that("Returned AUC from trainSLOPE is correct", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 100, q = 0.5, response = "binomial")
  x <- xy$x
  y <- xy$y

  tuned <- trainSLOPE(
    x,
    y,
    q = c(0.1, 0.2),
    measure = "auc",
    family = "binomial",
  )

  expect_equal(
    tuned$optima$mean[1],
    max(tuned$summary$mean)
  )
})
