test_that("plotting works", {
  set.seed(1)
  xy <- SLOPE:::randomProblem(100, 2)

  # one parameter
  fit <- SLOPE(xy$x, xy$y, alpha = 0.2)
  expect_silent(dont_plot(fit))

  # more parameters
  fit <- SLOPE(xy$x, xy$y, path_length = 10)
  expect_silent(dont_plot(fit))
})


test_that("plot.SLOPE works as expected", {
  fit <- SLOPE(heart$x, heart$y)
  p <- plot(fit)
  vdiffr::expect_doppelganger("plot.SLOPE-in-test", p)

  p <- plot(fit, intercept = TRUE, x_variable = "deviance_ratio")
  vdiffr::expect_doppelganger("plot.SLOPE-parameters-in-test", p)
})


test_that("plot.trainSLOPE works as expected", {
  set.seed(123)
  xy <- SLOPE:::randomProblem(1e2, 2)
  x <- xy$x
  y <- xy$y
  fit <- trainSLOPE(x, y, solver = "admm")

  expect_error(plot(fit, measure = "auc"))

  p <- plot(fit, measure = "mse")
  expect_s3_class(p, "ggplot")
  expect_silent(dont_plot(p))
  vdiffr::expect_doppelganger("trainSLOPE-in-test", p)

  fit <- trainSLOPE(subset(mtcars, select = c("mpg", "drat", "wt")),
                    mtcars$hp,
                    q = c(0.1, 0.2),
                    number = 10)
  p <- plot(fit, ci_col = "skyblue2", col = "black", ci_border = "mediumorchid")

  expect_s3_class(p, "ggplot")
  expect_silent(dont_plot(p))
  vdiffr::expect_doppelganger("trainSLOPE-double-q-in-test", p)

})
