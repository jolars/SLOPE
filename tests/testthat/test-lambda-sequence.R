test_that("lambda sequences are generated correctly", {
  for (type in c("gaussian", "bh", "oscar")) {
    # n = p
    lambda1 <- regularizationWeights(10, type = type, n = 10)

    # n > p
    lambda2 <- regularizationWeights(10, type = type, n = 100)

    # n < p
    lambda3 <- regularizationWeights(100, type = type, n = 10)

    expect_true(all(diff(lambda1) <= 0))
    expect_true(all(diff(lambda2) <= 0))
    expect_true(all(diff(lambda3) <= 0))
  }
})
