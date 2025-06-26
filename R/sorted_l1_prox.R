#' Sorted L1 Proximal Operator
#'
#' The proximal operator for the Sorted L1 Norm, which is the penalty function
#' in SLOPE. It solves the problem
#' \deqn{
#'   \arg\,\min_x
#'     \Big(J(x, \lambda) + \frac{1}{2} ||x - v||_2^2\Big)
#' }{
#'   argmin_x (J(x, \lambda) + 0.5||x - v||_2^2)
#' }
#' where \eqn{J(x, \lambda)} is the Sorted L1 Norm.
#'
#' @param x A vector. In SLOPE, this is the vector of coefficients.
#' @param lambda A non-negative and decreasing sequence
#'   of weights for the Sorted L1 Norm. Needs to be the same length as
#'   `x`.
#' @param method DEPRECATED
#' @return An evaluation of the proximal operator at `x` and `lambda`.
#'
#' @source
#' M. Bogdan, E. van den Berg, Chiara Sabatti, Weijie Su, and Emmanuel J.
#' Candès, “SLOPE – adaptive variable selection via convex optimization,” Ann
#' Appl Stat, vol. 9, no. 3, pp. 1103–1140, 2015.
#'
#' @export
sortedL1Prox <- function(x, lambda, method) {
  if (!missing(method)) {
    warning(
      "The 'method' argument is deprecated and ",
      "has no effect. It will be removed in a future version."
    )
  }

  stopifnot(
    length(x) == length(lambda),
    !is.unsorted(rev(lambda)),
    all(lambda >= 0),
    all(is.finite(lambda)),
    all(is.finite(x))
  )

  res <- sortedL1ProxCpp(as.matrix(x), lambda)

  as.vector(res)
}
