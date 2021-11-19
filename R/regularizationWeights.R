#' Generate Regularization (Penalty) Weights for SLOPE
#'
#' Please see [SLOPE()] for detailed information regarding the parameters in
#' this function.
#'
#' Note that these sequences are automatically scaled (unless a value for
#' the `alpha` parameter is manually supplied) when using [SLOPE()]. In this
#' function, nu such scaling is attempted.
#'
#' @param n_lambda The number of lambdas to generate. This should typically
#'   be equal to the number of predictors in your data set.
#' @param n The number of rows (observations) in the design matrix.
#' @inheritParams SLOPE
#'
#' @return A vector of length `n_lambda` with regularization weights.
#'
#' @seealso [SLOPE()]
#'
#' @export
#'
#' @examples
#' regularizationWeights(100, q = 0.2, lambda_type = "oscar")
regularizationWeights <- function(n_lambda,
                                  q = 0.2,
                                  theta1 = 1,
                                  theta2 = 0.5,
                                  lambda_type = c("bh", "gaussian", "oscar"),
                                  n = NULL) {
  stopifnot(
    is.numeric(n_lambda),
    is.numeric(q),
    is.numeric(theta1),
    is.numeric(theta2),
    n_lambda > 0,
    q >= 0,
    q <= 1,
    theta1 >= 0,
    theta2 >= 0
  )

  lambda_type <- match.arg(lambda_type)

  if (!is.null(n)) {
    stopifnot(is.numeric(n), floor(n) > 0)
  } else {
    n <- 1 # cpp function needs an integer value
  }

  if (lambda_type == "gaussian" && is.null(n)) {
    stop("`n` cannot be NULL when `lambda_type = 'gaussian''")
  }

  lambda <- lambdaSequence(n_lambda, q, theta1, theta2, lambda_type, n)

  as.vector(lambda)
}
