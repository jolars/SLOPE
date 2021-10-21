# test_that("model training works with trainSLOPE", {
#   set.seed(48)
#
#   for (family in c("gaussian", "binomial")) {
#     xy <- SLOPE:::randomProblem(1e2, 2, response = family)
#
#     fit <- trainSLOPE(xy$x,
#                       xy$y,
#                       family = family,
#                       n_folds = 2,
#                       q = c(0.1, 0.2),
#                       repeats = 2,
#                       path_length = 2)
#     expect_s3_class(fit, "TrainedSLOPE")
#   }
# })

test_that("plot.trainSLOPE works as expected", {
  xy <- SLOPE:::randomProblem(1e2, 2)
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, solver = "admm")

  expect_error(plot(fit, measure = "auc"))
  p <- plot(fit)
  expect_s3_class(p, "trellis")
  expect_silent(dont_plot(p))

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), n_folds = 2, solver = "admm")

  p <- plot(fit)
  expect_s3_class(p, "trellis")
})

test_that("Misclassification Rate works properly", {
  # set.seed(42)
  # xy <- SLOPE:::randomProblem(200, p=100, q=0.5, response="binomial")
  # x <- xy$x
  # y <- xy$y
  #
  # fit <- trainSLOPE(x, y, q = c(0.1, 0.2), n_folds = 2, measure = "misclass", family = "binomial")
  #
  # expect_equal(fit$measure$measure, "misclass")
  # expect_equal(fit$measure$label, "Misclassification Rate")
  # expect_equal(fit$optima$mean, 0.235)
})


test_that("trainSLOPE  returns error in the case of invalid measures", {

  set.seed(42)
  xy <- SLOPE:::randomProblem(200, p=100, q=0.5, response="gaussian")
  x <- xy$x
  y <- xy$y

  expect_error(trainSLOPE(x, y, q = c(0.1, 0.2), measure = "misclass", family = "gaussian"),
               "For the given family: gaussian, measure needs to be one of: mse, mae")

})










