test_that("penalty strength is invariant to number of observations", {
  x <- abalone$x
  y <- abalone$y

  x_small <- x
  y_small <- y
  x_large <- do.call(rbind, rep(list(x), 2))
  y_large <- rep(y, 2)

  for (family in c("gaussian", "poisson")) {
    for (scale in c("none", "sd", "l1", "l2")) {

      alpha <- switch(scale,
                      none = 0.1,
                      sd = 0.3,
                      l1 = 0.5,
                      l2 = 0.3)

      for (center in c(TRUE, FALSE)) {
        f0 <- SLOPE(x_small, y_small, family = family,
                    alpha = alpha, scale = scale, center = center,
                    tol_rel_gap = 1e-6, tol_infeas = 1e-6)
        f1 <- SLOPE(x_large, y_large, family = family,
                    alpha = alpha, scale = scale, center = center,
                    tol_rel_gap = 1e-6, tol_infeas = 1e-6)

        expect_equivalent(coef(f0), coef(f1), tol = 1e-3)
      }
    }
  }
})
