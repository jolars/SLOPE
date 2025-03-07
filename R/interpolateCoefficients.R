#' Interpolate coefficients
#'
#' @param interpolation_list a list generated from [interpolatePenalty()]
#' @param beta coefficients
#' @param intercepts intercepts
#'
#' @return A matrix (or list of matrices) with new coefficients based
#'   on linearly interpolating from new and old lambda values.
#' @keywords internal
interpolateCoefficients <- function(beta, intercepts, interpolation_list) {
  d <- length(interpolation_list$frac)

  ip_beta <- vector("list", d)
  ip_intercepts <- vector("list", d)

  for (i in seq_along(interpolation_list$left)) {
    left <- interpolation_list$left[i]
    right <- interpolation_list$right[i]
    frac <- interpolation_list$frac[i]
    ip_beta[[i]] <- beta[[left]] * frac + beta[[right]] * (1 - frac)
    ip_intercepts[[i]] <-
      intercepts[[left]] * frac + intercepts[[right]] * (1 - frac)
  }

  list(beta = ip_beta, intercepts = ip_intercepts)
}
