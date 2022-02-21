#' Plot results from diagnostics collected during model fitting
#'
#' This function plots various diagnostics collected during
#' the model fitting resulting from a call to [SLOPE()] *provided that
#' `diagnostics = TRUE`*.
#'
#' @importFrom stats reshape
#' @importFrom ggplot2 element_blank
#'
#' @param object an object of class `"SLOPE"`.
#' @param ind either "last"
#' @param xvar what to place on the x axis. `iteration` plots each iteration,
#'   `time` plots the wall-clock time.
#' @param ... parameters that will be used in \link[ggplot2]{theme}
#'
#' @return An object of class `"ggplot"`, which will be plotted on the
#'   current device unless stored in a variable.
#'
#' @seealso [SLOPE()], [ggplot2::theme()]
#'
#' @export
#'
#' @examples
#' x <- SLOPE(abalone$x, abalone$y, diagnostics = TRUE)
#' plotDiagnostics(x)
plotDiagnostics <- function(object,
                            ind = max(object$diagnostics$penalty),
                            xvar = c("time", "iteration"),
                            ...) {

  stopifnot(inherits(object, "SLOPE"),
            is.numeric(ind),
            length(ind) == 1)

  xvar <- match.arg(xvar)

  if (is.null(object$diagnostics))
    stop("no diagnostics found in fit;",
         "did you call SLOPE() with diagnostics = TRUE?")

  d <- object$diagnostics

  d <- subset(d, subset = d$penalty == ind)

  d <- reshape(d,
               direction = "long",
               varying = c("primal", "dual"),
               v.names = "Value",
               idvar = c("iteration", "time", "penalty"),
               timevar = "Variable",
               times = c("primal", "dual"))

  if (xvar == "time") {
    plt <- ggplot(d, aes(x = !!quote(time),
                         y = !!quote(Value),
                         col = !!quote(Variable))) +
      xlab("Time (seconds)")

  } else if (xvar == "iteration") {
    plt <- ggplot(d, aes(x = !!quote(iteration),
                         y = !!quote(Value),
                         col = !!quote(Variable))) +
      xlab("Iteration")
  }

  if (nrow(d) <= 1) {
    plt <- plt + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
  }

  plt <- plt +
    geom_line() +
    ylab("Objective") +
    theme(legend.title = element_blank(), ...)

  plt
}
