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
                      path_length = 2)
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
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p=100, q=0.5, response="binomial")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, family = "binomial")

  expect_equal(fit$measure$measure,
               c("mse", "mae", "deviance", "misclass", "auc"))
  expect_equal(fit$measure$label,
               c("Mean Squared Error", "Mean Absolute Error",
                 "Binomial Deviance", "Misclassification Rate", "AUC"))

  expect_equal(fit$optima, structure(list(q = c(0.1, 0.2, 0.2, 0.1, 0.2),
                                          alpha = c(0.0407419692328029,
                                                    0.0095161892230067,
                                                    7.46792437250563e-05,
                                                    0.00586052810915448,
                                                    0.0095161892230067),
                                          measure =
                                            c("auc", "deviance", "mae",
                                              "misclass", "mse"),
                                          mean = c(0.660558191808191,
                                                   0.980242385828498,
                                                   0.505232740237056,
                                                   0.235,
                                                   0.323146578355248),
                                          se = c(0.0427697302697303,
                                                 0.0895079147155515,
                                                 0.011858410680709,
                                                 0.055,
                                                 0.0303254835186373),
                                          lo = c(0.11711724249003,
                                                 -0.157063504055364,
                                                 0.354557346282327,
                                                 -0.463841260489608,
                                                 -0.0621752239560481),
                                          hi = c(1.20399914112635,
                                                 2.11754827571236,
                                                 0.655908134191785,
                                                 0.933841260489608,
                                                 0.708468380666544)),
                                     row.names = c(NA, -5L),
                                     class = "data.frame"),
               tolerance = 0.0001)
})


test_that("trainSLOPE works properly for gaussian family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p = 100, q = 0.5, response = "gaussian")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, family = "gaussian")

  expect_equal(fit$measure$measure,
               c("mse", "mae"))
  expect_equal(fit$measure$label,
               c("Mean Squared Error", "Mean Absolute Error"))

  expect_equal(fit$optima, structure(list(q = c(0.2, 0.2),
                                          alpha = c(0.060637891022855,
                                                    0.060637891022855),
                                          measure = c("mae", "mse"),
                                          mean = c(2.86335499483679,
                                                   13.1331062209287),
                                          se = c(0.172286938561511,
                                                 2.15038344415019),
                                          lo = c(0.67424188010548,
                                                 -14.1901060817242),
                                          hi = c(5.05246810956809,
                                                 40.4563185235815)),
                                     row.names = c(NA, -2L),
                                     class = "data.frame"),
               tolerance = 0.0001)
})

test_that("trainSLOPE works properly for poisson family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p=100, q=0.5, response="poisson")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, family = "poisson")

  expect_equal(fit$measure$measure, c("mse", "mae"))
  expect_equal(fit$measure$label, c("Mean Squared Error",
                                    "Mean Absolute Error"))

  expect_equal(fit$optima, structure(list(q = c(0.2, 0.1),
                                          alpha = c(2822597.89902365,
                                                    51737675.0334651),
                                          measure = c("mae", "mse"),
                                          mean = c(89508659.6408427,
                                                   803377308062305024),
                                          se = c(7450055.16497561,
                                                 131811729988185280),
                                          lo = c(-5153266.58113311,
                                                 -871449519796954624),
                                          hi = c(184170585.862818,
                                                 2478204135921564672)),
                                     row.names = c(NA, -2L),
                                     class = "data.frame"),
               tolerance = 0.0001)
})


test_that("trainSLOPE works properly for multinomial family", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p=100, q=0.5, response="multinomial")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, family = "multinomial")

  expect_equal(fit$measure$measure, c("mse", "mae", "deviance", "misclass"))
  expect_equal(fit$measure$label, c("Mean Squared Error",
                                    "Mean Absolute Error",
                                    "Multinomial Deviance",
                                    "Misclassification Rate"))

  expect_equal(fit$optima, structure(list(q = c(0.2, 0.2, 0.2, 0.2),
                                          alpha = c(0.0232836583559346,
                                                    6.93003501141371e-05,
                                                    0.000481771041719308,
                                                    0.0143391993455749),
                                          measure = c("deviance", "mae",
                                                      "misclass", "mse"),
                                          mean = c(212.520373211118,
                                                   0.323759418425716,
                                                   0.48,
                                                   0.216079211804543),
                                          se = c(13.5614867355816,
                                                 0.0312753402429364,
                                                 0.05,
                                                 0.0207765280569875),
                                          lo = c(40.2053462219011,
                                                 -0.0736314578945576,
                                                 -0.155310236808735,
                                                 -0.0479116073944179),
                                          hi = c(384.835400200335,
                                                 0.72115029474599,
                                                 1.11531023680873,
                                                 0.480070031003504)),
                                     row.names = c(NA, -4L),
                                     class = "data.frame"),
               tolerance = 0.0001)
})



test_that("trainSLOPE returns error in the case of invalid measures", {

  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p=100, q=0.5, response="gaussian")
  x <- xy$x
  y <- xy$y

  expect_error(trainSLOPE(x, y, q = c(0.1, 0.2),
                          measure = "misclass",
                          family = "gaussian"),
               "For the given family: gaussian,
               measure needs to be one of: mse, mae")
})










