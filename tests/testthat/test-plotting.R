test_that("plotting works", {
  set.seed(1)
  xy <- SLOPE:::randomProblem(100, 5)

  # one parameter
  fit <- SLOPE(xy$x, xy$y, alpha = 0.2)
  expect_silent(dont_plot(fit))

  # more parameters
  fit <- SLOPE(xy$x, xy$y, path_length = 100)
  expect_silent(dont_plot(fit))
})

test_that("plot.SLOPE works as expected", {
  fit <- SLOPE(heart$x, heart$y)
  expect_silent(dont_plot(fit))
  expect_silent(dont_plot(
    fit,
    intercept = TRUE,
    x_variable = "deviance_ratio"
  ))

  fit <- SLOPE(wine$x, wine$y, family = "multinomial")
  expect_silent(dont_plot(fit))
})

test_that("plot.trainedSLOPE works as expected", {
  set.seed(123)

  tune <- cvSLOPE(
    subset(mtcars, select = c("mpg", "drat", "wt")),
    mtcars$hp,
    q = c(0.1, 0.2),
    n_folds = 10
  )
  expect_silent(dont_plot(tune, ci_col = "salmon"))

  tune <- cvSLOPE(
    subset(mtcars, select = c("mpg", "drat", "wt")),
    mtcars$hp,
    q = 0.4,
    n_folds = 10
  )
  expect_silent(dont_plot(tune, ci_col = "salmon"))
})
