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

  alpha <- 0.01

  coef_ref <-
    c(
      -0.075761257437997,
      -0.0397870591523659,
      0.179401937344632,
      0.0719298289970407,
      -0.143050053903786,
      0.0119708915618286,
      -0.0152132541332644,
      0.0299541615757605,
      -0.0557280347240528,
      -0.0502826103600737,
      -0.0405067293862495
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

  expect_equivalent(coef_ref, as.matrix(coef(slope_fit)), tol = 1e-3)
})
