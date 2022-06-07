test_that("verbosity returns printouts", {
  xy <- SLOPE:::randomProblem(100, 10)
  x <- xy$x
  y <- xy$y

  expect_output(print(SLOPE(x, y, verbosity = 1)))
  expect_output(print(SLOPE(x, y, verbosity = 2)))
  expect_output(print(SLOPE(x, y, verbosity = 3)))
})

test_that("printing trainedSLOPE", {
  set.seed(1)
  tune <- trainSLOPE(subset(mtcars, select = c("mpg", "drat", "wt")),
    mtcars$hp,
    q = c(0.1, 0.2),
    number = 8,
    repeats = 5,
    measure = "mse"
  )
  expect_output(print(tune))
})


test_that("printing ABSLOPE", {
  set.seed(17)
  xy <- SLOPE:::randomProblem(100, 200, response = "gaussian")
  x <- as.matrix(xy$x)
  y <- xy$y
  fit <- ABSLOPE(x, y)
  expect_output(print(fit))
})
