# Compute the usual unbiased estimate of the variance in a linear model.
estimateNoise <- function(x, y, intercept = TRUE) {
  n <- NROW(x)
  p <- NCOL(x)

  stopifnot(n > p)

  # TODO: Don't do this
  if (intercept) {
    x <- cbind(1, x)
  }

  if (intercept) {
    fit <- stats::lm(y ~ x)
  } else {
    fit <- stats::lm(y ~ x - 1)
  }

  sqrt(sum(fit$residuals^2) / (n - p + intercept))
}
