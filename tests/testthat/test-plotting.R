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

test_that("plot.trainedSLOPE works as expected", {
  set.seed(123)
  tune <- trainSLOPE(
    subset(mtcars, select = c("mpg", "drat", "wt")),
    mtcars$hp,
    q = c(0.1, 0.2),
    number = 10
  )
  p <- plot(tune, ci_col = "salmon")
  vdiffr::expect_doppelganger("plot_trainedSLOPE-in-test", p)

  tune <- trainSLOPE(subset(mtcars, select = c("mpg", "drat", "wt")),
    mtcars$hp,
    q = 0.4,
    number = 10
  )

  p <- plot(tune, ci_col = "salmon")
  vdiffr::expect_doppelganger("q_plot_trainedSLOPE-in-test", p)

  xy <- SLOPE:::randomProblem(200, p = 10, q = 0.5, response = "binomial")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, family = "binomial")

  p <- plot(fit, ci_col = "salmon")

  vdiffr::expect_doppelganger("binom_plot_trainedSLOPE-in-test", p)
})
