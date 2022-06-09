#' Interpolate penalty values
#'
#' @param penalty current penalty sequence
#' @param x new sequence
#'
#' @return Interpolated values of lambda
#' @author Jerome Friedman, Trevor Hastie, Rob Tibshirani, and Noah Simon
#'
#' @keywords internal
interpolatePenalty <- function(penalty, x) {
  if (length(penalty) == 1) {
    nums <- length(x)
    left <- rep(1, nums)
    right <- left
    xfrac <- rep(1, nums)
  } else {
    x[x > max(penalty)] <- max(penalty)
    x[x < min(penalty)] <- min(penalty)
    k <- length(penalty)
    xfrac <- (penalty[1] - x) / (penalty[1] - penalty[k])
    penalty <- (penalty[1] - penalty) / (penalty[1] - penalty[k])
    coord <- stats::approx(penalty, seq(penalty), xfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    xfrac <- (xfrac - penalty[right]) / (penalty[left] - penalty[right])
    xfrac[left == right] <- 1
    xfrac[abs(penalty[left] - penalty[right]) < .Machine$double.eps] <- 1
  }

  list(
    left = left,
    right = right,
    frac = xfrac
  )
}
