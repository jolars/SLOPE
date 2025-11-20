#' Refit SLOPE Model with Optimal Parameters
#'
#' Refits a SLOPE model using the optimal parameters found through
#' cross-validation. This is a convenience function to avoid having to manually
#' extract optimal parameters and refit.
#'
#' @param object an object of class `'TrainedSLOPE'`, typically from a call to
#'   [cvSLOPE()] or [trainSLOPE()]
#' @param x the design matrix
#' @param y the response vector
#' @param measure which performance measure to use for selecting optimal
#'   parameters. If `NULL` (default), uses the first measure in the
#'   `TrainedSLOPE` object.
#' @param ... additional arguments passed to [SLOPE()]
#'
#' @return An object of class `'SLOPE'` fit with the optimal parameters
#'
#' @seealso [SLOPE()]
#' @family model-tuning
#'
#' @examples
#' # Cross-validation
#' tune <- trainSLOPE(
#'   bodyfat$x,
#'   bodyfat$y,
#'   q = c(0.1, 0.2),
#'   measure = "mse"
#' )
#'
#' # Refit with optimal parameters
#' fit <- refit(tune, bodyfat$x, bodyfat$y)
#'
#' # Use the fitted model
#' coef(fit)
#' predict(fit, bodyfat$x)
#'
#' @export
refit <- function(object, x, y, measure = NULL, ...) {
  UseMethod("refit")
}

#' @export
refit.TrainedSLOPE <- function(object, x, y, measure = NULL, ...) {
  if (is.null(measure)) {
    measure <- object$measure$measure[1]
  }

  if (!(measure %in% object$measure$measure)) {
    stop(
      "measure '",
      measure,
      "' not found in TrainedSLOPE object. ",
      "Available measures: ",
      paste(object$measure$measure, collapse = ", ")
    )
  }

  # Extract optimal parameters for the specified measure
  optima_row <- object$optima[object$optima$measure == measure, ]

  if (nrow(optima_row) == 0) {
    stop("No optimal parameters found for measure '", measure, "'")
  }

  q_opt <- optima_row$q[1]
  alpha_opt <- optima_row$alpha[1]

  # Extract gamma if present
  if ("gamma" %in% names(optima_row)) {
    gamma_opt <- optima_row$gamma[1]
  } else {
    gamma_opt <- 0
  }

  # Get additional arguments from the original call if available
  original_call <- as.list(object$call)

  # Build argument list for SLOPE
  slope_args <- list(
    x = x,
    y = y,
    q = q_opt,
    alpha = alpha_opt,
    gamma = gamma_opt
  )

  # Merge with user-provided arguments (user arguments take precedence)
  slope_args <- utils::modifyList(slope_args, list(...))

  # Call SLOPE with optimal parameters
  do.call(SLOPE, slope_args)
}
