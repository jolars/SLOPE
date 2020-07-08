test_that("different solvers return equivalent output", {
  set.seed(853)

  admm <- SLOPE(bodyfat$x, bodyfat$y, solver = "admm",
                tol_rel = 1e-6, tol_abs = 1e-7)
  fista <- SLOPE(bodyfat$x, bodyfat$y, solver = "fista",
                 tol_rel_gap = 1e-9)

  expect_equivalent(coef(admm), coef(fista), tol = 1e-3)
})
