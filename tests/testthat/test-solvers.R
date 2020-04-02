test_that("different solvers return equivalent output", {
  set.seed(853)

  admm <- SLOPE(bodyfat$x, bodyfat$y, solver = "admm")
  fista <- SLOPE(bodyfat$x, bodyfat$y, solver = "fista")

  expect_equivalent(coef(admm), coef(fista), tol = 1e-3)
})
