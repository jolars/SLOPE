test_that("different screening algorithms return equivalent results", {

  set.seed(1119)

  xy <- SLOPE:::randomProblem(50, 100)
  x <- xy$x
  y <- xy$y

  fit_strong <- SLOPE(x, y, screen_alg = "strong")
  fit_previous <- SLOPE(x, y, screen_alg = "previous")

  expect_equivalent(coef(fit_strong), coef(fit_previous), tol = 1e-3)
})
