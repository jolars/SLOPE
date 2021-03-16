test_that("incorrect data dimensions throw errors", {
  expect_error(SLOPE(matrix(1, 3, 3), double(2)))
})

test_that("na values in input throws errors", {
  expect_error(SLOPE(matrix(NA, 3, 3), double(3)))
  expect_error(SLOPE(matrix(3, 3, 3), c(NA, NA, 1)))
})

test_that("erroneous lambda ans alpha input throws errors", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(SLOPE(x, y, lambda = 1:2))
  expect_error(SLOPE(x, y, lambda = 1:3))
  expect_error(SLOPE(x, y, lambda = -c(1, 2, 3)))
  expect_error(SLOPE(x, y, alpha = -1))
  expect_error(SLOPE(x, y, alpha = 1:2))
})

test_that("erroneous standardization settings throw", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(SLOPE(x, y, scale = 3))
})

test_that("unsupported error settings are caught", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(SLOPE(x, y, family = "binomial", solver = "admm"))
})

test_that("sparse matrix and centering throws", {
  xy <- SLOPE:::randomProblem(density = 0.1)
  expect_error(SLOPE(xy$x, xy$y, center = TRUE))
})

test_that("specifying alpha estimation throws if not gaussian", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(SLOPE(x, y, alpha = "estimate", family = "binomial",
                     path_length = 1))
})

test_that("alpha estimation and path_length > 1 returns warning", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_warning(SLOPE(x, y, alpha = "estimate", path_length = 3))
})
