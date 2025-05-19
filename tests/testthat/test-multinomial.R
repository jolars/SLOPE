test_that("glmnet and SLOPE return same unpenalized model", {
  # Fixed predictors
  x1 <- c(
    1.2, -0.5, 0.8, -1.1, 0.3, 1.5, -0.2, 0.7, -0.9, 0.4,
    0.1, -1.3, 0.6, -0.7, 1.1, -0.4, 0.9, -1.0, 0.5, -0.8
  )
  x2 <- c(
    -0.3, 0.7, -1.2, 0.4, -0.8, 0.2, -0.5, 1.1, -0.9, 0.6,
    -1.0, 0.3, -0.7, 0.8, -0.4, 1.3, -0.6, 0.5, -1.1, 0.9
  )

  # Fixed response (deliberately creating a pattern)
  y <- c(
    1, 2, 3, 1, 2, 3, 1, 2, 3, 1,
    2, 3, 1, 2, 3, 1, 2, 3, 1, 2
  )

  x <- scale(cbind(x1, x2))
  y <- factor(y)

  g_coef <- structure(c(
    0.201106343886854, 0.286875567305304, 0.395516767951309,
    0.175083912785424, 0.200156047570187, 0.631090659835371
  ), dim = 3:2, dimnames = list(
    c("", "x1", "x2"), c("1", "2")
  ))

  ofit <- SLOPE(x, y, family = "multinomial", alpha = 1e-9)
  coefs <- as.matrix(do.call(cbind, coef(ofit)))

  expect_equivalent(g_coef, coefs, tol = 1e-4)
})
