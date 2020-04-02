test_that("deprecated arguments and functions return warnings", {
  xy <- SLOPE:::randomProblem(100, 10)
  x <- xy$x
  y <- xy$y

  expect_warning(SLOPE(x, y, solver = "matlab"))
  expect_warning(SLOPE(x, y, normalize = TRUE))
  expect_warning(SLOPE(x, y, fdr = 0.1))
  expect_warning(SLOPE(X = x, y = y))
})
