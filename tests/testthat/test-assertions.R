test_that("incorrect data dimensions throw errors", {
  expect_error(SLOPE(matrix(1, 3, 3), double(2)))
})

test_that("na values in input throws errors", {
  expect_error(SLOPE(matrix(NA, 3, 3), double(3)))
  expect_error(SLOPE(matrix(3, 3, 3), c(NA, NA, 1)))
})

test_that("erroneous lambda and alpha input throws errors", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(SLOPE(x, y, lambda = 1:2))
  expect_error(SLOPE(x, y, lambda = 1:3))
  expect_error(SLOPE(x, y, lambda = -c(1, 2, 3)))
  expect_error(SLOPE(x, y, alpha = -1))
})

test_that("erroneous standardization settings throw", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(SLOPE(x, y, scale = 3))
})

test_that("unsupported error settings are caught", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(SLOPE(x, y, family = "binomial"))
})

test_that("specifying alpha estimation throws if not gaussian", {
  x <- matrix(1, 3, 3)
  y <- double(3)

  expect_error(
    SLOPE(x, y, alpha = "estimate", family = "binomial", path_length = 1)
  )
})

test_that("SLOPE returns an error for empty `y`", {
  expect_error(SLOPE(x = NULL, y = NULL, q = 0.1), "`y` is empty")
})


test_that("SLOPE returns an error for `alpha='estimate'` and `family='gaussian'`", {
  expect_error(
    SLOPE(
      bodyfat$x,
      bodyfat$y,
      alpha = "estimate",
      family = "poisson"
    ),
    "`alpha = 'estimate'` can only be used if `family = 'gaussian'`"
  )
})
