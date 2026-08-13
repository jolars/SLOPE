test_that("patterns", {
  set.seed(8)

  x <- as.matrix(abalone$x)
  y <- abalone$y

  g <- SLOPE(x, y, family = "gaussian", alpha = 1e-4, patterns = TRUE)

  pattern <- g[["patterns"]]

  expect_equal(nrow(pattern[[1]]), ncol(x))
})
