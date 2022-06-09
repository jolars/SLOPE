#' Plot results from diagnostics collected during model fitting
#'
#' This function plots various diagnostics collected during
#' the model fitting resulting from a call to [SLOPE()] *provided that
#' `diagnostics = TRUE`*.
#'
#' @param object an object of class `"SLOPE"`.
#' @param ind either "last"
#' @param xvar what to place on the x axis. `iteration` plots each iteration,
#'   `time` plots the wall-clock time.
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
                            xvar = c("time", "iteration")) {
  stopifnot(
    inherits(object, "SLOPE"),
    is.numeric(ind),
    length(ind) == 1
  )

  xvar <- match.arg(xvar)

  if (is.null(object$diagnostics)) {
    stop(
      "no diagnostics found in fit; ",
      "did you call SLOPE() with diagnostics = TRUE?"
    )
  }

  d <- object$diagnostics

  d <- subset(d, subset = d$penalty == ind)

  d <- stats::reshape(
    d,
    direction = "long",
    varying = c("primal", "dual"),
    v.names = "Value",
    idvar = c("iteration", "time", "penalty"),
    timevar = "Variable",
    times = c("primal", "dual")
  )

  if (xvar == "time") {
    plt <- ggplot2::ggplot(
      d, ggplot2::aes(
        x = !!quote(time),
        y = !!quote(Value),
        col = !!quote(Variable)
      )
    ) +
      ggplot2::xlab("Time (seconds)")
  } else if (xvar == "iteration") {
    plt <- ggplot2::ggplot(
      d, ggplot2::aes(
        x = !!quote(iteration),
        y = !!quote(Value),
        col = !!quote(Variable)
      )
    ) +
      ggplot2::xlab("Iteration")
  }

  if (nrow(d) <= 1) {
    plt <- plt + ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  }

  plt <- plt +
    ggplot2::geom_line() +
    ggplot2::ylab("Objective") +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  plt
}
