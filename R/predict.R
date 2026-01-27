#' Generate Predictions from SLOPE Models
#'
#' Return predictions from models fit by [SLOPE()].
#'
#' @param object an object of class `"SLOPE"`, typically the result of
#'   a call to [SLOPE()]
#' @param x new data
#' @param type type of prediction; `"link"` returns the linear predictors,
#'   `"response"` returns the result of applying the link function,
#'    and `"class"` returns class predictions.
#' @param ... ignored and only here for method consistency
#' @param alpha penalty parameter for SLOPE models; if `NULL`, the
#'   values used in the original fit will be used
#' @param simplify if `TRUE`, [base::drop()] will be called before returning
#'   the coefficients to drop extraneous dimensions
#' @param exact if `TRUE` and the given parameter values differ from those in
#'   the original fit, the model will be refit by calling [stats::update()] on
#'   the object with the new parameters. If `FALSE`, the predicted values
#'   will be based on interpolated coefficients from the original
#'   penalty path.
#' @param sigma deprecated. Please use `alpha` instead.
#'
#' @seealso [stats::predict()], [stats::predict.glm()]
#' @family SLOPE-methods
#'
#' @return Predictions from the model with scale determined by `type`.
#'
#' @export
#'
#' @examples
#' fit <- with(mtcars, SLOPE(cbind(mpg, hp), vs, family = "binomial"))
#' predict(fit, with(mtcars, cbind(mpg, hp)), type = "class")
predict.SLOPE <- function(
  object,
  x,
  alpha = NULL,
  type = "link",
  simplify = TRUE,
  exact = FALSE,
  sigma,
  ...
) {
  # This method (the base method) only generates linear predictors

  if (inherits(x, "sparseMatrix")) {
    x <- as_dgCMatrix(x)
  }

  if (inherits(x, "data.frame")) {
    x <- as.matrix(x)
  }

  if (!missing(sigma)) {
    warning("`sigma` is deprecated. Please use `alpha` instead.")
    alpha <- sigma
  }

  # Get coefficients with intercepts using coef()
  has_intercept <- getElement(object, "has_intercept")
  coefs <- stats::coef(
    object,
    alpha = alpha,
    exact = exact,
    simplify = FALSE,
    intercept = has_intercept,
    ...
  )

  # Extract beta and intercepts from each matrix in the list
  n_penalties <- length(coefs)
  beta <- vector("list", n_penalties)
  intercepts <- vector("list", n_penalties)

  for (i in seq_len(n_penalties)) {
    if (has_intercept) {
      intercepts[[i]] <- coefs[[i]][1, ]
      beta[[i]] <- coefs[[i]][-1, , drop = FALSE]
    } else {
      intercepts[[i]] <- rep(0, NCOL(coefs[[i]]))
      beta[[i]] <- coefs[[i]]
    }
  }

  n <- NROW(x)
  p <- NROW(beta[[1L]])
  m <- NCOL(beta[[1L]])

  stopifnot(p == NCOL(x))

  lin_pred <- array(
    dim = c(n, m, n_penalties),
    dimnames = list(
      rownames(x),
      dimnames(beta[[1L]])[[2]]
    )
  )

  for (i in seq_len(n_penalties)) {
    lp <- as.matrix(x %*% beta[[i]])
    lin_pred[,, i] <- sweep(lp, 2, intercepts[[i]], "+") # nolint
  }

  lin_pred
}

#' @rdname predict.SLOPE
#' @export
predict.GaussianSLOPE <- function(
  object,
  x,
  sigma = NULL,
  type = c("link", "response"),
  simplify = TRUE,
  ...
) {
  type <- match.arg(type)

  out <- NextMethod(object, type = type) # always linear predictors

  if (simplify) {
    out <- drop(out)
  }

  out
}

#' @rdname predict.SLOPE
#' @export
predict.BinomialSLOPE <- function(
  object,
  x,
  sigma = NULL,
  type = c("link", "response", "class"),
  simplify = TRUE,
  ...
) {
  type <- match.arg(type)

  lin_pred <- NextMethod(object, type = type)

  out <- switch(
    type,
    link = lin_pred,
    response = 1 / (1 + exp(-lin_pred)),
    class = {
      cnum <- ifelse(lin_pred > 0, 2, 1)
      clet <- object$class_names[cnum]

      if (is.matrix(cnum)) {
        clet <- array(clet, dim(cnum), dimnames(cnum))
      }

      clet
    }
  )

  if (simplify) {
    out <- drop(out)
  }

  out
}

#' @rdname predict.SLOPE
#' @export
predict.PoissonSLOPE <- function(
  object,
  x,
  sigma = NULL,
  type = c("link", "response"),
  exact = FALSE,
  simplify = TRUE,
  ...
) {
  type <- match.arg(type)

  lin_pred <- NextMethod(object, type = type)

  out <- switch(type, link = lin_pred, response = exp(lin_pred))

  if (simplify) {
    out <- drop(out)
  }

  out
}

#' @export
#' @rdname predict.SLOPE
predict.MultinomialSLOPE <- function(
  object,
  x,
  sigma = NULL,
  type = c("link", "response", "class"),
  exact = FALSE,
  simplify = TRUE,
  ...
) {
  type <- match.arg(type)

  lin_pred <- NextMethod(object, type = type, simplify = FALSE)
  m <- NCOL(lin_pred)

  out <- switch(
    type,
    response = {
      n <- nrow(lin_pred)
      m <- ncol(lin_pred)
      path_length <- dim(lin_pred)[3]

      tmp <- array(0, c(n, m, path_length))
      tmp[, 1:m, ] <- lin_pred

      aperm(
        apply(
          tmp,
          c(1, 3),
          function(x) exp(x) / sum(1 + exp(x))
        ),
        c(2, 1, 3)
      )
    },
    link = lin_pred,
    class = {
      response <- stats::predict(object, x, type = "response", simplify = FALSE)
      tmp <- apply(response, c(1, 3), which.max)
      class_names <- object$class_names

      predicted_classes <-
        apply(
          tmp,
          2,
          function(a) factor(a, levels = 1:(m + 1), labels = class_names)
        )
      colnames(predicted_classes) <- dimnames(lin_pred)[[3]]
      predicted_classes
    }
  )

  if (simplify) {
    out <- drop(out)
  }

  out
}
