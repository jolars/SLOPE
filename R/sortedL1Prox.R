#' Sorted L1 Proximal Operator
#'
#' The proximal operator for the Sorted L1 Norm, which is the penalty function
#' in SLOPE. It solves the problem
#' \deqn{
#'   \arg\,\min_x \Big( J(x, \lambda) + \frac{1}{2} |x - v|_2^2\Big)
#' }{
#'   argmin_x (J(x, \lambda) + 0.5|x - v|_2^2)
#' }
#' where \eqn{J(x, \lambda)} is the Sorted L1 Norm.
#'
#' @param x A vector. In SLOPE, this is the vector of coefficients.
#' @param lambda A non-negative and decreasing sequence
#'   of weights for the Sorted L1 Norm. Needs to be the same length as
#'   `x`.
#' @param method Method used in the prox. `stack` is a stack-based algorithm
#'   (Algorithm 4 in Bogdan et al.).
#' @return The
#'
#' @source
#' M. Bogdan, E. van den Berg, Chiara Sabatti, Weijie Su, and Emmanuel J.
#' Candès, “SLOPE – adaptive variable selection via convex optimization,” Ann
#' Appl Stat, vol. 9, no. 3, pp. 1103–1140, 2015, doi: 10.1214/15-AOAS842.
#'
#' @export
sortedL1Prox <- function(x, lambda) {
  stopifnot(
    length(x) == length(lambda),
    !is.unsorted(rev(lambda)),
    all(lambda >= 0),
    all(is.finite(lambda)),
    all(is.finite(x))
  )

  res <- sorted_l1_prox(as.matrix(x), lambda)
  as.vector(res)
}


