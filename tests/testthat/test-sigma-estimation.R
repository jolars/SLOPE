test_that("sigma estimation works", {
  set.seed(808)

  xy <- SLOPE:::randomProblem(n = 1000, p = 10)

  expect_silent(SLOPE(xy$x, xy$y, sigma = "estimate", n_sigma = 1))

  # large p
  xy <- SLOPE:::randomProblem(n = 100, p = 250)

  expect_silent(SLOPE(xy$x, xy$y, sigma = "estimate", n_sigma = 1))
})
