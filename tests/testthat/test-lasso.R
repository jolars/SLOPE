test_that("lasso and slope fits are equivalent if all lambda are equal", {
  set.seed(1)
  xy <- SLOPE:::randomProblem(100, 10)
  x <- xy$x
  y <- xy$y

  lambda <- 0.6

  lasso_coef <-
    new(
      "dgCMatrix",
      i = c(0L, 2L, 3L),
      p = c(0L, 3L),
      Dim = c(
        11L,
        1L
      ),
      Dimnames = list(
        c(
          "(Intercept)",
          "V1",
          "V2",
          "V3",
          "V4",
          "V5",
          "V6",
          "V7",
          "V8",
          "V9",
          "V10"
        ),
        "s0"
      ),
      x = c(
        -0.0560909669856042,
        2.25670882779904,
        -2.41538332234389
      ),
      factors = list()
    )

  slope <- SLOPE(
    x,
    y,
    center = FALSE,
    scale = FALSE,
    lambda = rep(lambda, ncol(x)),
    alpha = 1
  )

  slope_coef <- coef(slope)

  expect_equivalent(lasso_coef, slope_coef, tol = 1e-3)
})
