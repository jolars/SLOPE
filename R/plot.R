#' Plot coefficients
#'
#' Plot the fitted model's regression
#' coefficients along the regularization path.
#'
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap xlab ylab theme labs
#'
#' @param x an object of class `"SLOPE"`
#' @param ... parameters that will be used in \link[ggplot2]{theme}
#' @param intercept whether to plot the intercept
#' @param x_variable what to plot on the x axis. `"alpha"` plots
#'   the scaling parameter for the sequence, `"deviance_ratio"` plots
#'   the fraction of deviance explained, and `"step"` plots step number.
#'
#' @seealso [ggplot2::theme()], [SLOPE()], [plotDiagnostics()]
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
                       x_variable = c("alpha", "deviance_ratio", "step"),
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
              step = seq_along(object[["alpha"]]))

  xlab <- switch(x_variable,
                 alpha = expression(alpha),
                 deviance_ratio = "Fraction Deviance Explained",
                 step = "Step")


  d <- as.data.frame(as.table(coefs))
  d[["x"]] <- rep(x, each = p*m)

  if(m > 1) {
    ggplot(d, aes(x = !!quote(x),
                  y = !!quote(Freq),
                  col = !!quote(Var1))) +
      geom_line() +
      facet_wrap(~!!quote(Var2)) +
      ylab(expression(hat(beta))) +
      xlab(xlab) +
      labs(color = 'Variable name') +
      theme(...)
  }else {
    ggplot(d, aes(x = !!quote(x), y = !!quote(Freq), col = !!quote(Var1))) +
      geom_line() +
      ylab(expression(hat(beta))) +
      xlab(xlab) +
      labs(color = 'Variable name') +
      theme(...)
  }
}

#' Plot results from cross-validation
#'
#' @importFrom ggplot2 geom_vline geom_ribbon
#' @importFrom grDevices colors
#'
#' @param x an object of class `'TrainedSLOPE'`, typically from a call
#'   to [trainSLOPE()]
#' @param measure any of the measures used in the call to [trainSLOPE()]. If
#'   `measure = "auto"` then deviance will be used for binomial and multinomial
#'   models, whilst mean-squared error will be used for Gaussian and Poisson
#'   models.
#' @param ci_alpha alpha (opacity) for fill in confidence limits
#' @param ci_col color for border of confidence limits
#' @param ... other arguments that are passed on to [lattice::xyplot()]
#' @param plot_min whether to mark the location of the penalty corresponding
#'   to the best prediction score
#' @param ci_border color (or flag to turn off and on) the border of the
#'   confidence limits
#'
#' @seealso [trainSLOPE()], [lattice::xyplot()], [lattice::panel.xyplot()]
#' @family model-tuning
#'
#' @return An object of class `"trellis"` is returned and, if used
#'   interactively, will most likely have its print function
#'   [lattice::print.trellis()] invoked, which draws the plot on the
#'   current display device.
#'
#' @export
#'
#' @examples
#' # Cross-validation for a SLOPE binomial model
#' set.seed(123)
#' tune <- trainSLOPE(subset(mtcars, select = c("mpg", "drat", "wt")),
#'                    mtcars$hp,
#'                    q = c(0.1, 0.2),
#'                    number = 10)
#' plot(tune, ci_col = "salmon", col = "black")
plot.TrainedSLOPE <- function(x,
                              measure = c("auto", "mse", "mae",
                                          "deviance", "auc", "misclass"),
                              plot_min = TRUE,
                              ci_alpha = 0.2,
                              ci_border = FALSE,
                              ci_col = "salmon",
                              ...) {

  if(!(ci_col %in% colors()))
    stop("ci_col", ci_col, "is not a valid color representation.")

  if(!(ci_border %in% colors() | is.logical(ci_border)))
    stop("ci_border is", ci_border, "when it should be logical or a valid
           color representation.")

  object <- x
  family <- object[["model"]][["family"]]

  measure <- match.arg(measure)

  if (measure == "auto") {
    measure <- switch(
      family,
      gaussian = "mse",
      binomial = "deviance",
      multinomial = "deviance",
      poisson = "mse"
    )

    if (!any(measure %in% object[["measure"]][["measure"]]))
      measure <- object[["measure"]][["measure"]][1]
  }

  if (!any(measure %in% object[["measure"]][["measure"]]))
    stop("measure ", measure, " was not used or not available when",
         "fitting the model")

  if (length(measure) > 1)
    stop("you are only allowed to plot one measure at a time")


  measure_label <- object[["measure"]][["label"]][object[["measure"]][["measure"]] == measure]

  summary <- object[["summary"]][object[["summary"]][["measure"]] == measure,]
  optimum <- object[["optima"]][object[["optima"]][["measure"]] == measure, ,
                                drop = FALSE]
  optimum[["label_q"]] <- paste0("q = ", as.factor(optimum[["q"]]))
  model <- object[["model"]]

  q <- unique(summary[["q"]])

  summary[["q"]] <- as.factor(summary[["q"]])
  summary[["label_q"]] <- paste0("q = ", as.factor(summary[["q"]]))

  xlab <- expression(log[e](alpha))

  p <- ggplot(summary, aes(x = log(!!quote(alpha)), y = mean)) +
    geom_line() +
    xlab(xlab) +
    ylab(measure_label)

  if (length(q) > 1) {

    p <- p + facet_wrap(~label_q)
  }

  if(plot_min) {

    p <- p + geom_vline(data = optimum,
                        aes(xintercept = log(!!quote(alpha))),
                        linetype = "dotted")
  }

  border_col <- c(ci_border, NA, ci_col)[is.character(ci_border) +
                                           2*isFALSE(ci_border) +
                                           3*isTRUE(ci_border)]

  p <- p + geom_ribbon(aes(ymin = !!quote(lo), ymax = !!quote(hi)),
                       fill = ci_col,
                       color = border_col,
                       alpha = ci_alpha)

  p
}
