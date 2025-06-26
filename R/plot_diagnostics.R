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
#' @return Invisibly returns NULL. The function is called for its
#'   side effect of producing a plot.
#'
#' @seealso [SLOPE()]
#'
#' @export
#'
#' @examples
#' x <- SLOPE(abalone$x, abalone$y, diagnostics = TRUE)
#' plotDiagnostics(x)
plotDiagnostics <- function(
  object,
  ind = max(object$diagnostics$penalty),
  xvar = c("time", "iteration")
) {
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

  primal <- d[["primal"]]
  dual <- d[["dual"]]

  primal_color <- "#0072B2" # Blue
  dual_color <- "#D55E00" # Vermillion

  # Set up the plot window
  if (xvar == "time") {
    x_primal <- d$time
    x_dual <- d$time
    xlab <- "Time (seconds)"
  } else if (xvar == "iteration") {
    x_primal <- d$iteration
    x_dual <- d$iteration
    xlab <- "Iteration"
  }

  # Create the plot
  plot(
    x_primal,
    primal,
    type = "n",
    xlab = xlab,
    ylab = "Objective",
    ylim = range(c(primal, dual))
  )

  # Add lines for primal and dual
  graphics::lines(x_primal, primal, col = primal_color)
  graphics::lines(x_dual, dual, col = dual_color)

  # Add legend
  if (nrow(d) > 1) {
    graphics::legend(
      "topright",
      legend = c("primal", "dual"),
      col = c(primal_color, dual_color)
    )
  }

  invisible()
}
