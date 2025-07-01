test_that("deviance() works.", {
  fit <- SLOPE(abalone$x, abalone$y, family = "poisson", path_length = 20)
  d <- deviance(fit)
  coef_ref <- c(
    0.0243668140665241,
    -0.0634487290701952,
    0,
    0.250962832580353,
    1.29093061706281,
    0,
    0,
    0,
    0.493869779302109
  )
  expect_equivalent(
    as.vector(fit$coefficients[[5]]),
    coef_ref,
    tolerance = 1e-3
  )
})
