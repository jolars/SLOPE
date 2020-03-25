test_that("sigma estimation works", {
  xy <- SLOPE:::randomProblem(n = 1000, p = 10)

  expect_silent(SLOPE(xy$x, xy$y, sigma = "estimate", n_sigma = 1))
})
