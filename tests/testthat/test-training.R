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

test_that("plot.trainSLOPE works as expected", {
  xy <- SLOPE:::randomProblem(1e2, 2)
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, solver = "admm")

  expect_error(plot(fit, measure = "auc"))
  p <- plot(fit)
  expect_s3_class(p, "trellis")
  expect_silent(dont_plot(p))

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, solver = "admm")

  p <- plot(fit)
  expect_s3_class(p, "trellis")
})


test_that("trainSLOPE works properly for binomial family", {
  set.seed(1454)
  xy <- SLOPE:::randomProblem(200, p = 10, q = 0.5, response = "binomial")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, family = "binomial")

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
    fit$optima,
    structure(list(
      q = c(0.1, 0.1, 0.2, 0.1, 0.1),
      alpha = c(
        0.0857618386463366,
        0.00177453086249557,
        8.57618386463366e-06,
        0.0200315767754891,
        0.00177453086249557
      ),
      measure = c("auc", "deviance", "mae", "misclass", "mse"),
      mean = c(
        0.82236720470667,
        0.357468195162962,
        0.179541901444967,
        0.07,
        0.111500400822224
      ),
      se = c(
        0.0395490593835067,
        0.0591266351433017,
        0.0364692338159319,
        0.01,
        0.017110631194972
      ),
      lo = c(
        0.319848759056704,
        -0.393806936328931,
        -0.283843649991689,
        -0.057062047361747,
        -0.105910782306268
      ),
      hi = c(
        1.32488565035664,
        1.10874332665486,
        0.642927452881622,
        0.197062047361747,
        0.328911583950716
      )
    ), row.names = c(NA, -5L), class = "data.frame"),
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
    tol_rel_gap = 1e-5,
    tol_infeas = 1e-7
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
    fit$optima,
    structure(list(
      q = c(0.2, 0.2),
      alpha = c(0.0291535803838849, 0.0291535803838849),
      measure = c("mae", "mse"),
      mean = c(0.815910923266303, 1.04527630688046),
      se = c(0.0275194184023733, 0.0566349151161424),
      lo = c(0.466243558825295, 0.325661480198884),
      hi = c(1.16557828770731, 1.76489113356204)
    ), row.names = c(NA, -2L), class = "data.frame"),
    tolerance = 0.001
  )
})

test_that("trainSLOPE works properly for poisson family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 5, q = 0.5, response = "poisson")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, family = "poisson")

  expect_equal(fit$measure$measure, c("mse", "mae"))
  expect_equal(
    fit$measure$label,
    c("Mean Squared Error", "Mean Absolute Error")
  )

  expect_equal(
    fit$optima,
    structure(list(
      q = c(0.2, 0.2),
      alpha = c(0.694874382791202, 0.162303301420616),
      measure = c("mae", "mse"),
      mean = c(2.17509028617307, 31.7907389590011),
      se = c(0.443936208583228, 20.1886058597683),
      lo = c(-3.46565406988658, -224.729820433151),
      hi = c(7.81583464223272, 288.311298351153)
    ), row.names = c(NA, -2L), class = "data.frame"),
    tolerance = 0.0001
  )
})


test_that("trainSLOPE works properly for multinomial family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 5, q = 0.5, response = "multinomial")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, family = "multinomial")

  expect_equal(fit$measure$measure, c("mse", "mae", "deviance", "misclass"))
  expect_equal(fit$measure$label, c(
    "Mean Squared Error",
    "Mean Absolute Error",
    "Multinomial Deviance",
    "Misclassification Rate"
  ))

  expect_equal(
    fit$optima,
    structure(list(
      q = c(0.1, 0.2, 0.1, 0.1),
      alpha = c(
        0.000271276009840165,
        9.11438098491086e-06,
        0.00497243693348613,
        3.90216946049708e-05
      ),
      measure = c("deviance", "mae", "misclass", "mse"),
      mean = c(55.1549471708103, 0.103076394047732, 0.125, 0.0548364722790797),
      se = c(1.20617967890997, 0.0063949508132034, 0.005, 0.000395127349159509),
      lo = c(
        39.8289812219668,
        0.0218208397374032,
        0.0614689763191265,
        0.0498159032837969
      ),
      hi = c(
        70.4809131196538,
        0.184331948358062,
        0.188531023680874,
        0.0598570412743624
      )
    ), row.names = c(NA, -4L), class = "data.frame"),
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
      family = "gaussian"
    ),
    "For the given family: gaussian, measure needs to be one of: mse, mae"
  )
})
