test_that("incorrect data dimensions throw errors", {
  expect_error(SLOPE(matrix(1, 3, 3), double(2)))
})

test_that("na values in input throws errors", {
  expect_error(SLOPE(matrix(NA, 3, 3), double(3)))
  expect_error(SLOPE(matrix(3, 3, 3), c(NA, NA, 1)))
})

test_that("erroneous lambda ans sigma input throws errors", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(SLOPE(x, y, lambda = 1:2))
  expect_error(SLOPE(x, y, lambda = 1:3))
  expect_error(SLOPE(x, y, lambda = -c(1, 2, 3)))
  expect_error(SLOPE(x, y, sigma = -1))
  expect_error(SLOPE(x, y, sigma = 1:2))
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

test_that("specifying sigma estimation throws if not gaussian", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(SLOPE(x, y, sigma = "estimate", family = "binomial",
                     n_sigma = 1))
})

test_that("sigma estimation and n_sigma > 1 returns warning", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_warning(SLOPE(x, y, sigma = "estimate", n_sigma = 3))
})
