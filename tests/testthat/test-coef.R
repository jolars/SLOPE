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

test_that("coefficient rescaling", {
  set.seed(1150)
  xy <- SLOPE:::randomProblem(100, 10)

  xy$x[, 1] <- xy$x[, 1] * 10

  # check for slope
  fit <- SLOPE(xy$x, xy$y, alpha = 0.5)

  # check simplify
  coefs_original <- coef(fit, scale = "original")
  coefs_unscaled <- coef(fit, scale = "normalized")

  expect_false(all(coefs_original[-1, ] == coefs_unscaled[-1, ]))

  sds <- apply(xy$x, 2, function(x) {
    sqrt((length(x) - 1) / length(x)) * stats::sd(x)
  })

  expect_equal(
    as.vector(coefs_original[-1, ] * sds),
    as.vector(coefs_unscaled[-1, ]),
    tolerance = 1e-6
  )
})
