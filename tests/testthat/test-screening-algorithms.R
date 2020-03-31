test_that("different screening algorithms return equivalent results", {

  set.seed(1119)

  xy <- SLOPE:::randomProblem()
  x <- xy$x
  y <- xy$y

  fit_strong <- SLOPE(x, y, screen_alg = "strong")
  fit_working <- SLOPE(x, y, screen_alg = "working")

  expect_equivalent(coef(fit_strong), coef(fit_working), tol = 1e-3)
})
