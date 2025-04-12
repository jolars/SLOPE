test_that("patterns", {

  set.seed(8)

  x <- as.matrix(abalone$x)
  y <- abalone$y

  lm_fit <- lm(y ~ as.matrix(x))

  g <- SLOPE(x, y, family = "gaussian", alpha = 1e-4, patterns = TRUE)

  pattern <- g[["patterns"]]

  expect_equal(nrow(pattern[[1]]), ncol(x))
})
