test_that("unregularized logistic regression matches output from glm()", {
  set.seed(1)
  x <- scale(matrix(rnorm(3000), ncol = 3))
  x1 <- x[, 1]
  x2 <- x[, 2]
  x3 <- x[, 3]
  z <- 1 + 2 * x1 + 3 * x2 + x3
  pr <- 1 / (1 + exp(-z))
  y <- rbinom(1000, 1, pr)

  df <- data.frame(y = y, x1 = x1, x2 = x2)
  glm_fit <- glm(y ~ x1 + x2 + x3, data = df, family = "binomial")

  g_model <- SLOPE(
    cbind(x1, x2, x3),
    y,
    family = "binomial",
    diagnostics = TRUE,
    alpha = 1e-7
  )

  intercept <- g_model$intercept

  expect_equivalent(
    coef(glm_fit),
    as.vector(coef(g_model)),
    tol = 1e-3
  )
})

test_that("regularized slope logistic regression picks out correct features", {
  set.seed(2)
  p <- 10
  n <- 200
  k <- 3

  x <- matrix(rnorm(p * n), n, p)
  beta <- double(p)
  nz <- sample(p, k)

  beta[nz] <- 10
  z <- x %*% beta + 1
  prob <- 1 / (1 + exp(-z))

  y <- rbinom(n, 1, prob)

  slope_fit <- SLOPE(x, y, family = "binomial", alpha = 1 / sqrt(n))

  expect_setequal(nz, which(slope_fit$nonzeros[[1]]))
})

test_that("logistic regression works on the glioma dataset", {
  set.seed(825)

  glioma_fit <- SLOPE(
    x = glioma$x,
    y = glioma$y,
    family = "binomial",
    alpha = 1e-2
  )

  y_hat <- predict(glioma_fit, glioma$x, type = "class")

  accuracy <- mean(y_hat == glioma$y)

  expect_equal(accuracy, 0.975, tolerance = 1e-1)
})
