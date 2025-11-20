#' Model Deviance
#'
#' @param object an object of class `'SLOPE'`.
#' @param ... ignored
#'
#' @return For Gaussian models this is twice the residual sums of squares. For
#'   all other models, two times the negative loglikelihood is returned.
#' @export
#'
#' @seealso [SLOPE()]
#' @family SLOPE-methods
#'
#' @examples
#' fit <- SLOPE(heart$x, heart$y, family = "binomial")
#' deviance(fit)
deviance.SLOPE <- function(object, ...) {
  deviance_ratio <- object$deviance_ratio
  null_deviance <- object$null_deviance

  (1 - deviance_ratio) * null_deviance
}
