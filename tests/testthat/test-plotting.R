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

test_that("single solution dot chart works for gaussian", {
  set.seed(123)
  xy <- SLOPE:::randomProblem(50, 10)
  
  # Single solution with default parameters
  fit <- SLOPE(xy$x, xy$y, alpha = 0.1)
  expect_silent(dont_plot(fit))
  
  # Single solution with intercept
  expect_silent(dont_plot(fit, intercept = TRUE))
  
  # Single solution with magnitudes
  expect_silent(dont_plot(fit, magnitudes = TRUE))
  
  # Single solution without mark_zero
  expect_silent(dont_plot(fit, mark_zero = FALSE))
})

test_that("single solution dot chart works for multinomial", {
  set.seed(456)
  
  # Single solution multinomial
  fit <- SLOPE(wine$x, wine$y, family = "multinomial", alpha = 0.05)
  expect_silent(dont_plot(fit))
  
  # With intercept
  expect_silent(dont_plot(fit, intercept = TRUE))
  
  # With magnitudes
  expect_silent(dont_plot(fit, magnitudes = TRUE))
})

test_that("single solution uses variable names when available", {
  set.seed(789)
  X <- matrix(rnorm(100), 20, 5)
  colnames(X) <- paste0("Gene_", 1:5)
  y <- rnorm(20)
  
  fit <- SLOPE(X, y, alpha = 0.1)
  expect_silent(dont_plot(fit))
})

test_that("single solution uses class names for multinomial", {
  # Use wine data which has class names
  fit <- SLOPE(wine$x, wine$y, family = "multinomial", alpha = 0.05)
  expect_silent(dont_plot(fit))
})

test_that("dotchart arguments are passed through correctly", {
  set.seed(321)
  xy <- SLOPE:::randomProblem(30, 8)
  fit <- SLOPE(xy$x, xy$y, alpha = 0.15)
  
  # Test custom arguments
  expect_silent(dont_plot(fit, pch = 16))
  expect_silent(dont_plot(fit, col = "red"))
  expect_silent(dont_plot(fit, xlab = "Custom Label"))
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
