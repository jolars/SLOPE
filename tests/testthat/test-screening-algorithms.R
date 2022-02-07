test_that("different screening algorithms return equivalent results", {

  set.seed(1119)

  xy <- SLOPE:::randomProblem(50, 100)
  x <- xy$x
  y <- xy$y

  fit_strong <- SLOPE(x, y, screen_alg = "strong", tol_dev_change = 0,
                      tol_dev_ratio = 2, max_variables = 500)
  fit_previous <- SLOPE(x, y, screen_alg = "previous", tol_dev_change = 0,
                        tol_dev_ratio = 2, max_variables = 500)

  expect_equivalent(coef(fit_strong), coef(fit_previous), tol = 6e-3)
})
