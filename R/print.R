#' Print Results from SLOPE Fit
#'
#' @param x an object of class `'SLOPE'` or `'TrainedSLOPE'`
#' @param ... other arguments passed to [print()]
#'
#' @return Prints output on the screen
#'
#' @examples
#' fit <- SLOPE(wine$x, wine$y, family = "multinomial")
#' print(fit, digits = 1)
#' @method print SLOPE
#' @family SLOPE-methods
#' @seealso [SLOPE()], [print.SLOPE()]
#'
#' @export
print.SLOPE <- function(x, ...) {
  alpha <- x$alpha
  n_nonzero <- vapply(x$nonzeros, sum, FUN.VALUE = double(1))
  deviance_ratio <- x$deviance_ratio

  out <- data.frame(
    alpha = alpha,
    deviance_ratio = deviance_ratio,
    n_nonzero = n_nonzero
  )

  # print call
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  # print path summary
  cat("Path summary:\n")
  print(out, ...)
}

#' @rdname print.SLOPE
#' @method print TrainedSLOPE
#' @export
print.TrainedSLOPE <- function(x, ...) {
  cat(
    "\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )

  cat("Optimum values:\n")

  print(x$optima, ...)
}
