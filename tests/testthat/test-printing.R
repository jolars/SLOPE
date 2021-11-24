test_that("verbosity returns printouts", {
  xy <- SLOPE:::randomProblem(100, 10)
  x <- xy$x
  y <- xy$y

  expect_output(print(SLOPE(x, y, verbosity = 1)))
  expect_output(print(SLOPE(x, y, verbosity = 2)))
  expect_output(print(SLOPE(x, y, verbosity = 3)))
})
