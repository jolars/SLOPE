
#' Plot coefficients
#'
#' Plot the fitted model's regression
#' coefficients along the regularization path.
#'
#' @param x an object of class `"SLOPE"`
#' @param intercept whether to plot the intercept
#' @param x_variable what to plot on the x axis. `"alpha"` plots
#'   the scaling parameter for the sequence, `"deviance_ratio"` plots
#'   the fraction of deviance explained, and `"step"` plots step number.
#' @param ... further arguments passed to or from other methods.
#'
#' @seealso [SLOPE()], [plotDiagnostics()]
#' @family SLOPE-methods
#'
#' @return An object of class `"ggplot"`, which will be plotted on the
#'   current device unless stored in a variable.
#' @export
#'
#' @examples
#' fit <- SLOPE(heart$x, heart$y)
#' plot(fit)
plot.SLOPE <- function(x,
                       intercept = FALSE,
                       x_variable = c(
                         "alpha",
                         "deviance_ratio",
                         "step"
                       ),
                       ...) {
  object <- x

  x_variable <- match.arg(x_variable)

  coefs <- getElement(object, "coefficients")

  intercept_in_model <- "(Intercept)" %in% rownames(coefs)
  include_intercept <- intercept && intercept_in_model

  if (include_intercept) {
    coefs <- coefs[, , , drop = FALSE]
  } else if (intercept_in_model) {
    coefs <- coefs[-1, , , drop = FALSE]
  }

  p <- NROW(coefs) # number of features
  m <- NCOL(coefs) # number of responses

  x <- switch(x_variable,
    alpha = object[["alpha"]],
    deviance_ratio = object[["deviance_ratio"]],
    step = seq_along(object[["alpha"]])
  )

  xlab <- switch(x_variable,
    alpha = expression(alpha),
    deviance_ratio = "Fraction Deviance Explained",
    step = "Step"
  )

  d <- as.data.frame(as.table(coefs))
  d[["x"]] <- rep(x, each = p * m)

  plt <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = !!quote(x),
      y = !!quote(Freq),
      col = !!quote(Var1)
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = xlab,
      y = expression(hat(beta)),
      color = NULL
    )

  if (x_variable == "alpha") {
    plt <- plt + ggplot2::scale_x_log10()
  }

  if (m > 1) {
    plt <- plt + ggplot2::facet_wrap("Var2")
  }

  plt
}

#' Plot results from cross-validation
#'
#' @param x an object of class `'TrainedSLOPE'`, typically from a call
#'   to [trainSLOPE()]
#' @param measure any of the measures used in the call to [trainSLOPE()]. If
#'   `measure = "auto"` then deviance will be used for binomial and multinomial
#'   models, whilst mean-squared error will be used for Gaussian and Poisson
#'   models.
#' @param ci_alpha alpha (opacity) for fill in confidence limits
#' @param ci_col color for border of confidence limits
#' @param plot_min whether to mark the location of the penalty corresponding
#'   to the best prediction score
#' @param ci_border color (or flag to turn off and on) the border of the
#'   confidence limits
#' @param ... words
#'
#' @seealso [trainSLOPE()]
#' @family model-tuning
#'
#' @return An object of class `"ggplot"`, which will be plotted on the
#'   current device unless stored in a variable.
#'
#' @export
#'
#' @examples
#' # Cross-validation for a SLOPE binomial model
#' set.seed(123)
#' tune <- trainSLOPE(subset(mtcars, select = c("mpg", "drat", "wt")),
#'   mtcars$hp,
#'   q = c(0.1, 0.2),
#'   number = 10
#' )
#' plot(tune, ci_col = "salmon")
plot.TrainedSLOPE <- function(x,
                              measure = c(
                                "auto",
                                "mse",
                                "mae",
                                "deviance",
                                "auc",
                                "misclass"
                              ),
                              plot_min = TRUE,
                              ci_alpha = 0.2,
                              ci_border = FALSE,
                              ci_col = "salmon",
                              ...) {
  if (!(ci_col %in% grDevices::colors())) {
    stop("ci_col", ci_col, "is not a valid color representation.")
  }

  if (!(ci_border %in% grDevices::colors() | is.logical(ci_border))) {
    stop(
      "ci_border is",
      ci_border,
      "when it should be logical or a validcolor representation."
    )
  }

  object <- x
  family <- object[["model"]][["family"]]

  measure <- match.arg(measure)

  if (measure == "auto") {
    measure <- switch(family,
      gaussian = "mse",
      binomial = "deviance",
      multinomial = "deviance",
      poisson = "mse"
    )

    if (!any(measure %in% object[["measure"]][["measure"]])) {
      measure <- object[["measure"]][["measure"]][1]
    }
  }

  if (!any(measure %in% object[["measure"]][["measure"]])) {
    stop(
      "measure ",
      measure,
      " was not used or not available when",
      "fitting the model"
    )
  }

  if (length(measure) > 1) {
    stop("you are only allowed to plot one measure at a time")
  }

  measure_label <-
    object[["measure"]][["label"]][object[["measure"]][["measure"]] == measure]
  summary <-
    object[["summary"]][object[["summary"]][["measure"]] == measure, ]
  opt <- object[["optima"]]
  optimum <- opt[opt[["measure"]] == measure, , drop = FALSE]

  q <- unique(summary[["q"]])

  xlab <- expression(alpha)

  border_col <-
    c(
      ci_border,
      NA,
      ci_col
    )[is.character(ci_border) + 2 * isFALSE(ci_border) + 3 * isTRUE(ci_border)]

  p <- ggplot2::ggplot(
    summary,
    ggplot2::aes(x = !!quote(alpha), y = mean)
  ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(measure_label) +
    ggplot2::scale_x_log10() +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = !!quote(lo), ymax = !!quote(hi)),
      fill = ci_col,
      color = border_col,
      alpha = ci_alpha
    ) +
    ggplot2::geom_line()

  if (length(q) > 1) {
    p <- p + ggplot2::facet_wrap(
      "q",
      labeller = ggplot2::labeller(q = ggplot2::label_both)
    )
  }

  if (plot_min) {
    p <- p + ggplot2::geom_vline(
      data = optimum,
      ggplot2::aes(xintercept = !!quote(alpha)),
      linetype = "dotted"
    )
  }

  p
}
