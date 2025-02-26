test_that("interpolating coefficients works properly", {
  set.seed(3)
  xy <- SLOPE:::randomProblem(100, 10)

  # check for slope
  fit <- SLOPE(xy$x, xy$y)

  expect_s4_class(coef(fit), "dgCMatrix")
  expect_silent(coef(fit, alpha = c(0.001, 0.04)))
  # coef(fit, alpha = c(0.001, 0.04))

  # check for lasso
  fit <- SLOPE(xy$x, xy$y)

  expect_silent(coef(fit))
  expect_silent(coef(fit, lambda = c(0.2, 20)))

  # penalties are in the path already
  expect_silent(coef(fit, lambda = fit$lambda[c(2, 3)]))
})

test_that("simplify argument in coef() works as expected", {
  set.seed(1623)
  xy <- SLOPE:::randomProblem(100, 10)

  # check for slope
  fit <- SLOPE(xy$x, xy$y)

  # check simplify
  coefs <- coef(fit, simplify = TRUE)
  expect_s4_class(coefs, "dgCMatrix")

  coefs <- coef(fit, simplify = FALSE)
  expect_type(coefs, "list")
})

test_that("refitting works if exact = TRUE", {
  set.seed(1624)
  xy <- SLOPE:::randomProblem(100, 10)

  # check for slope
  fit <- SLOPE(xy$x, xy$y)

  # check simplify
  coefs <- coef(fit, alpha = 0.4, exact = TRUE, x = xy$x, y = xy$y)

  expect_s4_class(coefs, "dgCMatrix")
})
