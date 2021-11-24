test_that("screening rules return correct results for instances with known violations", {
  set.seed(216) # there is a violation for this seed and this setup

  for (solver in c("fista", "admm")) {
    for (family in c("gaussian", "binomial", "poisson", "multinomial")) {

      if (family != "gaussian" && solver == "admm")
        next

      d <- SLOPE:::randomProblem(100, 10, q = 0.1, response = family)

      fit0 <- SLOPE(d$x, d$y, family = family, screen = FALSE)
      fit1 <- SLOPE(d$x, d$y, family = family, screen = TRUE)

      expect_equivalent(coef(fit0), coef(fit1), 1e-4)
    }
  }
})

test_that("basic screening rule works", {
  set.seed(213)

  xy <- SLOPE:::randomProblem(100, 10, q = 0.1)

  fit <- SLOPE(xy$x, xy$y, screen = TRUE)

  expect_lt(length(fit$active_sets[[1]]), ncol(xy$x))
})
