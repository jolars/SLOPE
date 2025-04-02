test_that("alpha estimation works", {
  set.seed(808)

  xy <- SLOPE:::randomProblem(n = 1000, p = 10)

  expect_silent(SLOPE(xy$x, xy$y, alpha = "estimate"))

  # large p
  xy <- SLOPE:::randomProblem(n = 100, p = 80)

  expect_silent(SLOPE(xy$x, xy$y, alpha = "estimate"))
})
