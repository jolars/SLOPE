#' Generate Regularization (Penalty) Weights for SLOPE
#'
#' This function generates sequences of regularizations weights for use in
#' [SLOPE()] (or elsewhere).
#'
#' Please see [SLOPE()] for detailed information regarding the parameters in
#' this function, in particular the section *Regularization Sequences*.
#'
#' Note that these sequences are automatically scaled (unless a value for
#' the `alpha` parameter is manually supplied) when using [SLOPE()]. In this
#' function, nu such scaling is attempted.
#'
#' @param n_lambda The number of lambdas to generate. This should typically
#'   be equal to the number of predictors in your data set.
#' @param n The number of rows (observations) in the design matrix.
#' @param type The type of lambda sequence to use. See documentation for
#'   in [SLOPE()], including that related to the `lambda` parameter in that
#'   function.
#' @inheritParams SLOPE
#'
#' @return A vector of length `n_lambda` with regularization weights.
#'
#' @seealso [SLOPE()]
#'
#' @export
#'
#' @examples
#' # compute different penalization sequences
#' bh <- regularizationWeights(100, q = 0.2, type = "bh")
#'
#' gaussian <- regularizationWeights(
#'   100,
#'   q = 0.2,
#'   n = 300,
#'   type = "gaussian"
#' )
#'
#' oscar <- regularizationWeights(
#'   100,
#'   theta1 = 1.284,
#'   theta2 = 0.0182,
#'   type = "oscar"
#' )
#'
#' lasso <- regularizationWeights(100, type = "lasso") * mean(oscar)
#'
#' # Plot a comparison between these sequences
#' plot(bh, type = "l", ylab = expression(lambda))
#' lines(gaussian, col = "dark orange")
#' lines(oscar, col = "navy")
#' lines(lasso, col = "red3")
#'
#' legend(
#'   "topright",
#'   legend = c("BH", "Gaussian", "OSCAR", "lasso"),
#'   col = c("black", "dark orange", "navy", "red3"),
#'   lty = 1
#' )
regularizationWeights <- function(n_lambda = 100,
                                  type = c("bh", "gaussian", "oscar", "lasso"),
                                  q = 0.2,
                                  theta1 = 1,
                                  theta2 = 0.5,
                                  n = NULL) {
  stopifnot(
    is.numeric(n_lambda),
    is.numeric(q),
    is.numeric(theta1),
    is.numeric(theta2),
    n_lambda > 0,
    q > 1e-6,
    q < 1,
    theta1 >= 0,
    theta2 >= 0,
    isTRUE(is.finite(theta1)),
    isTRUE(is.finite(theta2))
  )

  type <- match.arg(type)

  if (!is.null(n) && type == "gaussian") {
    stopifnot(is.numeric(n), is.finite(n), floor(n) > 0)
  } else {
    n <- 1 # cpp function needs an integer value
  }

  if (type == "gaussian" && is.null(n)) {
    stop("`n` cannot be NULL when `type = 'gaussian''")
  }

  lambda <- lambdaSequence(n_lambda, q, theta1, theta2, type, n)

  as.vector(lambda)
}
