# Compute the usual unbiased estimate of the variance in a linear model.
estimateNoise <- function(x, y, intercept = TRUE) {
  n <- NROW(x)
  p <- NCOL(x)

  stopifnot(n > p)

  fit <- stats::lm.fit(x, y)
  sqrt(sum(fit$residuals^2) / (n - p + intercept))
}
