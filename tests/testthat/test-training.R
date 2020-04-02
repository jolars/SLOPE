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
                      n_sigma = 2)
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
