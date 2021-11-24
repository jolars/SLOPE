test_that("Regularization weights are reasonable", {
  is_valid_sequence <- function(x) {
    !is.unsorted(rev(x)) && all(is.finite(x)) && all(x >= 0)
  }

  # BH sequence type
  for (q in c(1e-5, 0.1, 0.99)) {
    for (n_lambda in c(1, 2, 1000)) {
      x <- regularizationWeights(n_lambda, q = q, type = "bh")
      expect_true(is_valid_sequence(x))
    }
  }

  # Gaussian sequence type
  for (q in c(1e-5, 0.1, 0.99)) {
    for (n_lambda in c(1, 2, 1000)) {
      for (n in c(50, 2000)) {
        x <- regularizationWeights(n_lambda, q = q, type = "gaussian", n = n)
        expect_true(is_valid_sequence(x))
      }
    }
  }

  # OSCAR
  for (theta1 in c(0, 0.5, 1)) {
    for (theta2 in c(0, 0.5, 1)) {
      x <- regularizationWeights(
        n_lambda,
        theta1 = theta1,
        theta2 = theta2,
        type = "oscar"
      )
      expect_true(is_valid_sequence(x))
    }
  }

  # lasso
  for (n_lambda in c(1, 100)) {
    x <- regularizationWeights(n_lambda, type = "lasso")
    expect_true(is_valid_sequence(x))
  }
})
