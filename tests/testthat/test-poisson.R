test_that("output from unregularized poisson model matches glm", {
  set.seed(654)

  for (intercept in c(TRUE, FALSE)) {
    xy <- SLOPE:::randomProblem(100, 10, response = "poisson", amplitude = 1)
    x <- xy$x
    y <- xy$y

    x <- scale(x)

    slope_fit <- SLOPE(
      x,
      y,
      family = "poisson",
      alpha = 1e-5,
      intercept = intercept,
      center = FALSE,
      scale = "none"
    )

    glm_fit <- if (intercept) {
      glm(y ~ 1 + ., data = data.frame(y = y, x), family = "poisson")
    } else {
      glm(y ~ 0 + ., data = data.frame(y = y, x), family = "poisson")
    }

    expect_equivalent(coef(glm_fit), as.matrix(coef(slope_fit)), tol = 1e-4)
  }
})

test_that("SLOPE reproduces lasso fit when all lambda are equal", {
  set.seed(0978213)

  xy <- SLOPE:::randomProblem(100, 10, response = "poisson", amplitude = 1)
  x <- xy$x
  y <- xy$y

  x <- scale(x)

  n <- nrow(x)
  p <- ncol(x)

  alpha <- 1

  gnt_coef <-
    c(
      0.490003407560466,
      0,
      0.628217573067706,
      0,
      0,
      0,
      -0.562605853351947,
      0,
      0,
      0,
      0
    )

  slope_fit <- SLOPE(
    x,
    y,
    family = "poisson",
    alpha = alpha,
    lambda = rep(1, p),
    scale = "none",
    center = FALSE,
    tol = 1e-6
  )

  expect_equivalent(gnt_coef, as.matrix(coef(slope_fit)), tol = 1e-3)
})
