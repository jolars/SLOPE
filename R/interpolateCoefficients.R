#' Interpolate coefficients
#'
#' @param interpolation_list a list generated from [interpolatePenalty()]
#' @param beta coefficients
#'
#' @return A matrix (or list of matrices) with new coefficients based
#'   on linearly interpolating from new and old lambda values.
#' @keywords internal
interpolateCoefficients <- function(beta,
                                    interpolation_list) {

  d <- length(interpolation_list$frac)

  ip_beta <- array(NA,
                   c(dim(beta)[1], dim(beta)[2], d),
                   dimnames = list(rownames(beta),
                                   colnames(beta),
                                   paste(seq_len(d))))

  for (i in seq_along(interpolation_list$left)) {
    ip_beta[, , i] <-
      beta[, , interpolation_list$left[i]] * interpolation_list$frac[i] +
      beta[, , interpolation_list$right[i]] * (1 - interpolation_list$frac[i])
  }

  ip_beta
}
