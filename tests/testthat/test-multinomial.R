test_that("glmnet and SLOPE return same unpenalized model", {

  set.seed(129)

  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  prob <- matrix(c(rep(1, n),
                   exp(3 + 2*x1 + x2),
                   exp(-1 + x1 - 3*x2)),
                 ncol = 3)
  prob <- sweep(prob, 1, apply(prob, 1, sum), "/")

  y = double(n)

  for (i in 1:n)
    y[i] <- sample(3, 1, replace = TRUE, prob = prob[i, ])

  x <- scale(cbind(x1, x2))
  y <- factor(y)

  library(glmnet)
  fit <- glmnet(x, y, family = "multinomial", lambda = 0, thresh = 1e-10,
                standardize = FALSE)
  g_coef <- as.matrix(do.call(cbind, coef(fit)))

  g_coef[,] <- g_coef[,] - g_coef[, 3]
  g_coef <- g_coef[, 1:2]

  ofit <- SLOPE(x, y, family = "multinomial", alpha = 1e-9)

  expect_equivalent(g_coef, coef(ofit), tol = 1e-4)
})

test_that("trainSLOPE for multinomial case works properly", {
  set.seed(42)
  xy <- SLOPE:::randomProblem(100, p = 20, response = "multinomial")
  x <- xy$x
  y <- xy$y

  fit <- trainSLOPE(x, y, q = c(0.1, 0.2), number = 2, family = "multinomial")

  # expect_equal(fit$measure$measure, c("mse", "mae", "deviance"))
  # expect_equal(fit$optima$mean, c(0.032337127203516, 0.0714196560684897, 0.0940813752794205))
})


