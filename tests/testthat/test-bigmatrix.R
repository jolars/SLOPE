library(SLOPE)
library(testthat)
library(bigmemory)

test_that("Bigmatrix", {
  n <- 4

  x <- diag(n)
  y <- c(8, 6, 4, 2)
  lambda <- c(4, 3, 2, 1)

  tmp_dir <- tempdir()

  x <- bigmemory::filebacked.big.matrix(
    n,
    n,
    type = "double",
    init = 0,
    backingfile = "example.bin",
    descriptorfile = "example.desc",
    backingpath = tmp_dir
  )

  diag(x) <- 1

  res <- SLOPE(
    x,
    y,
    family = "gaussian",
    intercept = FALSE,
    center = FALSE,
    scale = "none",
    lambda = lambda / n,
    alpha = 1
  )

  unlink(tmp_dir, recursive = TRUE)

  beta <- coef(res)

  expect_equal(as.vector(beta), c(4, 3, 2, 1), check.attributes = FALSE)
})
