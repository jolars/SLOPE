test_that("penalty strength is invariant to number of observations", {
  x <- abalone$x
  y <- abalone$y

  x_small <- x
  y_small <- y
  x_large <- do.call(rbind, rep(list(x), 2))
  y_large <- rep(y, 2)

  for (family in c("gaussian", "poisson")) {
    for (scale in c("none", "sd", "l1", "l2")) {
      alpha <- switch(scale, none = 0.1, sd = 0.3, l1 = 0.5, l2 = 0.3)

      for (center in c(TRUE, FALSE)) {
        f0 <- SLOPE(
          x_small,
          y_small,
          family = family,
          alpha = alpha,
          scale = scale,
          center = center
        )
        f1 <- SLOPE(
          x_large,
          y_large,
          family = family,
          alpha = alpha,
          scale = scale,
          center = center
        )

        expect_equivalent(coef(f0), coef(f1), tol = 1e-3)
      }
    }
  }
})

test_that("relative alpha uses the weighted null-model penalty", {
  x <- cbind(1:10, sin(1:10))
  y <- 2 * x[, 1] - x[, 2]
  weights <- c(0, 1:9)
  max_alpha <- SLOPE(x, y, path_length = 2, weights = weights)$alpha[1]
  fractions <- c(1.1, 0.4, 0.1)

  fit <- do.call(
    SLOPE,
    list(
      x = x,
      y = y,
      alpha = fractions,
      weights = weights,
      alpha_scale = "relative"
    )
  )
  absolute <- SLOPE(x, y, alpha = fractions * max_alpha, weights = weights)

  expect_equal(fit$alpha, fractions * max_alpha)
  expect_equal(coef(fit), coef(absolute), tolerance = 1e-5)

  new_alpha <- 0.3 * max_alpha
  expected <- SLOPE(x, y, alpha = new_alpha, weights = weights)
  expect_equal(
    coef(fit, alpha = new_alpha, exact = TRUE),
    coef(expected),
    tolerance = 1e-5
  )
})
