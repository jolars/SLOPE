#' Obtain coefficients
#'
#' This function returns coefficients from a model fit by [SLOPE()].
#'
#' If `exact = FALSE` and `alpha` is not in `object`,
#' then the returned coefficients will be approximated by linear interpolation.
#' If coefficients from another type of penalty sequence
#' (with a different `lambda`) are required, however,
#' please use [SLOPE()] to refit the model.
#'
#' @param object an object of class `'SLOPE'`.
#' @param intercept whether to include the intercept in the output; only
#'   applicable when `simplify = TRUE` and an intercept has been fit.
#' @param scale whether to return the coefficients in the original scale
#'   or in the normalized scale.
#' @param ... arguments that are passed on to [stats::update()] (and therefore
#'   also to [SLOPE()]) if `exact = TRUE` and the given penalty
#'   is not in `object`
#' @inheritParams predict.SLOPE
#'
#' @seealso [predict.SLOPE()], [SLOPE()]
#' @family SLOPE-methods
#'
#' @return Coefficients from the model.
#'
#' @export
#' @examples
#' fit <- SLOPE(mtcars$mpg, mtcars$vs, path_length = 10)
#' coef(fit)
#' coef(fit, scale = "normalized")
coef.SLOPE <- function(
  object,
  alpha = NULL,
  exact = FALSE,
  simplify = TRUE,
  intercept = TRUE,
  scale = c("original", "normalized"),
  sigma,
  ...
) {
  if (!missing(sigma)) {
    warning("`sigma` is deprecated. Please use `alpha` instead.")
    alpha <- sigma
  }

  scale <- match.arg(scale)

  beta <- switch(
    scale,
    original = getElement(object, "coefficients"),
    normalized = getElement(object, "coefficients_scaled")
  )
  intercepts <- switch(
    scale,
    original = getElement(object, "intercepts"),
    normalized = getElement(object, "intercepts_scaled")
  )
  penalty <- object$alpha
  value <- alpha

  if (is.null(value)) {
  } else if (all(value %in% penalty)) {
    beta <- beta[[penalty %in% value]]
  } else if (exact) {
    object <- stats::update(object, alpha = alpha, ...)
    beta <- switch(
      scale,
      original = getElement(object, "coefficients"),
      normalized = getElement(object, "coefficients_scaled")
    )
    intercepts <- switch(
      scale,
      original = getElement(object, "intercepts"),
      normalized = getElement(object, "intercepts_scaled")
    )
  } else {
    stopifnot(value >= 0)
    interpolation_list <- interpolatePenalty(penalty, value)
    res <- interpolateCoefficients(beta, intercepts, interpolation_list)
    beta <- res$beta
    intercepts <- res$intercepts
  }

  m <- NCOL(beta[[1]])

  if (simplify) {
    has_intercept <- getElement(object, "has_intercept")

    beta_out <- vector("list", m)

    for (i in seq_len(m)) {
      beta_out[[i]] <- do.call(
        cbind,
        lapply(beta, function(x) x[, i, drop = FALSE])
      )

      if (intercept && has_intercept) {
        intercepts_i <- as.vector(do.call(cbind, lapply(intercepts, "[", i)))
        beta_out[[i]] <- rbind(
          intercepts_i,
          beta_out[[i]],
          deparse.level = 0
        )
      }
    }

    beta <- beta_out

    if (m == 1) {
      beta <- beta[[1]]
    }
  }

  beta
}
