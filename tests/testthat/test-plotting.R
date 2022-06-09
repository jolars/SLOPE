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
  p1 <- plot(fit)
  p2 <- plot(fit, intercept = TRUE, x_variable = "deviance_ratio")

  fit <- SLOPE(wine$x, wine$y, family = "multinomial")
  p3 <- plot(fit)

  skip_on_ci()

  vdiffr::expect_doppelganger("plot.SLOPE-in-test", p1)
  vdiffr::expect_doppelganger("plot.SLOPE-parameters-in-test", p2)
  vdiffr::expect_doppelganger("plot.SLOPE-multinomial-in-test", p3)
})

test_that("plot.trainedSLOPE works as expected", {
  set.seed(123)

  tune <- trainSLOPE(
    subset(mtcars, select = c("mpg", "drat", "wt")),
    mtcars$hp,
    q = c(0.1, 0.2),
    number = 10
  )
  p1 <- plot(tune, ci_col = "salmon")

  tune <- trainSLOPE(subset(mtcars, select = c("mpg", "drat", "wt")),
    mtcars$hp,
    q = 0.4,
    number = 10
  )
  p2 <- plot(tune, ci_col = "salmon")

  xy <- SLOPE:::randomProblem(200, p = 10, q = 0.5, response = "binomial")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, family = "binomial")

  p3 <- plot(fit, ci_col = "salmon")

  skip_on_ci()

  vdiffr::expect_doppelganger("plot_trainedSLOPE-in-test", p1)
  vdiffr::expect_doppelganger("q_plot_trainedSLOPE-in-test", p2)
  vdiffr::expect_doppelganger("binom_plot_trainedSLOPE-in-test", p3)
})
