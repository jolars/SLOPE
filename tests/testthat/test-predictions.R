test_that("predictions work for all models", {
  set.seed(1)

  for (family in c("gaussian", "binomial", "poisson", "multinomial")) {
    xy <- SLOPE:::randomProblem(100, 10, response = family)
    x <- xy$x
    y <- xy$y

    fit <- SLOPE(x, y, family = family, n_sigma = 5)

    for (type in c("link", "response", "class")) {
      if (type == "class" && family %in% c("gaussian", "poisson"))
        next

      expect_silent(predict(fit, x, type = type))
    }
  }
})

