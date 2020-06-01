test_that("penalty strength is invariant to number of observations", {
  set.seed(2026)

  xy <- SLOPE:::randomProblem(100, 3, q = 1)

  x <- xy$x
  y <- xy$y

  x <- abalone$x
  y <- abalone$y

  x_small <- do.call(rbind, rep(list(x), 10))
  y_small <- rep(y, 10)
  x_large <- x
  y_large <- y

  for (scale in c("none", "sd", "l1", "l2")) {

    alpha <- switch(scale,
                    none = 1.3,
                    sd = 0.3,
                    l1 = 1.5,
                    l2 = 1.2)

    for (center in c(TRUE, FALSE)) {
      f0 <- SLOPE(x_small, y_small,
                  alpha = alpha, scale = scale, center = center,
                  tol_rel_gap = 1e-8, tol_infeas = 1e-8)
      f1 <- SLOPE(x_large, y_large,
                  alpha = alpha, scale = scale, center = center,
                  tol_rel_gap = 1e-8, tol_infeas = 1e-8)

      expect_equivalent(coef(f0), coef(f1), tol = 1e-4)
    }
  }
})
