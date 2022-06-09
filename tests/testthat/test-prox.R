# Equivalent formulation for the SLOPE prox, but
# using isotonic regression (PAVA) algorithm instead
prox_sorted_l1_isotone <- function(x, lambda) {
  sign <- sign(x)
  x <- abs(x)
  ord <- order(x)

  res <- stats::isoreg(x[ord] - rev(lambda))

  out <- pmax.int(res$yf, 0)
  out[ord] <- out

  out * sign
}

test_that("Prox and isotonic regression agree", {
  library(SLOPE)
  n <- 15

  set.seed(2254)

  x <- rnorm(n)
  lambda <- sort(runif(n), decreasing = TRUE)

  out_isotone_ref <- prox_sorted_l1_isotone(x, lambda)
  out_stack <- sortedL1Prox(x, lambda, "stack")
  out_pava <- sortedL1Prox(x, lambda, "pava")

  tol <- .Machine$double.eps^0.95

  expect_equal(out_isotone_ref, out_stack)
  expect_equal(out_isotone_ref, out_pava)
})

test_that("Prox works for simple examples", {
  r1 <- sortedL1Prox(c(5, 2), c(4, 2))
  expect_identical(r1, c(1, 0))

  r2 <- sortedL1Prox(c(3, 3), c(0, 0))
  expect_identical(r2, c(3, 3))

  r3 <- sortedL1Prox(c(2, 1), c(3, 0))
  expect_identical(r3, c(0, 0))
})

test_that("Assertions for prox function", {
  expect_error(sortedL1Prox(1:3, 1:3))
  expect_error(sortedL1Prox(1:3, c(-1, 2, 1)))
  expect_error(sortedL1Prox(c(NA, 1, 2), 3:1))
  expect_error(sortedL1Prox(c(1, 1, 2), c(Inf, 3, 1)))
})
