test_that("lambda sequences are generated correctly", {
  for (type in c("gaussian", "bh", "oscar")) {
    # Casee n = p
    lambda1 <- regularizationWeights(10, type = type, n = 10)

    # Case n > p
    lambda2 <- regularizationWeights(10, type = type, n = 100)

    # Case n < p
    lambda3 <- regularizationWeights(100, type = type, n = 10)

    expect_true(all(diff(lambda1) <= 0))
    expect_true(all(diff(lambda2) <= 0))
    expect_true(all(diff(lambda3) <= 0))
  }
})
