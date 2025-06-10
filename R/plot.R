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
#' @param magnitudes whether to plot the magnitudes of the coefficients
#' @param ... arguments passed to [graphics::matplot()]
#'
#' @seealso [SLOPE()], [plotDiagnostics()]
#' @family SLOPE-methods
#'
#' @return Invisibly returns NULL. The function is called for its
#'   side effect of producing a plot.
#' @export
#'
#' @examples
#' fit <- SLOPE(heart$x, heart$y)
#' plot(fit)
plot.SLOPE <- function(
  x,
  intercept = FALSE,
  x_variable = c(
    "alpha",
    "deviance_ratio",
    "step"
  ),
  magnitudes = FALSE,
  ...
) {
  object <- x

  x_variable <- match.arg(x_variable)

  coefs <- getElement(object, "coefficients")

  include_intercept <- intercept && object$has_intercept

  m <- NCOL(coefs[[1L]]) # number of responses

  x <- switch(
    x_variable,
    alpha = object[["alpha"]],
    deviance_ratio = object[["deviance_ratio"]],
    step = seq_along(object[["alpha"]])
  )

  xlab <- switch(
    x_variable,
    alpha = expression(alpha),
    deviance_ratio = "Fraction Deviance Explained",
    step = "Step"
  )

  coefs <- stats::coef(object, simplify = TRUE, intercept = include_intercept)

  if (magnitudes) {
    coefs <- abs(coefs)
  }

  xlim <- if (x_variable == "alpha") rev(range(x)) else range(x)
  log_var <- if (x_variable == "alpha") "x" else ""

  default_plot_args <- list(
    type = "l",
    lty = 1,
    xlab = xlab,
    xlim = xlim,
    ylab = expression(hat(beta)),
    log = log_var,
    col = palette.colors(10, "Tableau")
  )

  plot_args <- utils::modifyList(
    default_plot_args,
    list(...)
  )

  if (m == 1) {
    xlim <- if (x_variable == "alpha") rev(range(x)) else range(x)
    plot_args_k <- utils::modifyList(plot_args, list(x = x, y = t(coefs)))
    do.call(graphics::matplot, plot_args_k)
  } else {
    for (k in seq_len(m)) {
      plot_args_k <- utils::modifyList(
        plot_args,
        list(
          x = x,
          y = t(coefs[[k]]),
          main = paste0("Response: ", k)
        )
      )
      do.call(graphics::matplot, plot_args_k)
    }
  }

  invisible()
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
#' @param plot_args list of additional arguments to pass to [plot()],
#'   which sets up the plot frame
#' @param polygon_args list of additional arguments to pass to
#'   [graphics::polygon()], which fills the confidence limits
#' @param lines_args list of additional arguments to pass to
#'  [graphics::lines()], which plots the mean
#' @param abline_args list of additional arguments to pass to
#'   [graphics::abline()], which plots the minimum
#' @param ... ignored
#'
#' @seealso [trainSLOPE()]
#' @family model-tuning
#'
#' @return A plot for every value of `q` is produced on the current device.
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
plot.TrainedSLOPE <- function(
  x,
  plot_min = TRUE,
  ci_alpha = 0.2,
  ci_border = NA,
  ci_col = "salmon",
  plot_args = list(),
  polygon_args = list(),
  lines_args = list(),
  abline_args = list(),
  measure,
  ...
) {
  object <- x
  family <- object[["model"]][["family"]]

  if (!missing(measure)) {
    warning(
      "`measure` is deprecated, and will be removed in a future version. ",
      "the measure will instead be taken from the `TrainedSLOPE` object",
      call. = FALSE
    )
  }

  measure <- object[["measure"]][["measure"]][1]
  measure_label <- object[["measure"]][["label"]][1]
  summary <- object[["summary"]]
  opt <- object[["optima"]]
  optimum <- opt[opt[["measure"]] == measure, , drop = FALSE]

  q <- unique(summary[["q"]])
  gamma <- unique(summary[["gamma"]])

  xlab <- expression(alpha)

  ci_col <- grDevices::adjustcolor(ci_col, alpha = ci_alpha)

  ylim <- range(c(summary[["lo"]], summary[["hi"]]))

  plot_args <- utils::modifyList(
    list(
      xlab = xlab,
      ylab = measure_label,
      log = "x",
      type = "n",
      ylim = ylim
    ),
    plot_args
  )

  polygon_args <- utils::modifyList(
    list(
      col = ci_col,
      border = ci_border
    ),
    polygon_args
  )

  lines_args <- utils::modifyList(list(), lines_args)

  abline_args <- utils::modifyList(
    list(
      v = optimum[["alpha"]],
      lty = 2
    ),
    abline_args
  )

  grid <- expand.grid(
    q = q,
    gamma = gamma
  )

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]

    ind <- which(
      summary[["q"]] == g[["q"]] &
        summary[["gamma"]] == g[["gamma"]]
    )
    summary_i <- summary[ind, ]

    # Plot frame
    plot_args <- utils::modifyList(
      plot_args,
      list(
        x = summary_i[["alpha"]],
        y = summary_i[["mean"]]
      )
    )

    do.call(plot, plot_args)

    if (nrow(grid) > 1) {
      title_parts <- NULL
      if (length(q) > 1) {
        q_part <- bquote(paste("q = ", .(g[["q"]])))
        title_parts <- if (is.null(title_parts)) q_part else
          bquote(paste(.(title_parts), ", ", .(q_part)))
      }
      if (length(gamma) > 1) {
        gamma_part <- bquote(paste(gamma, " = ", .(g[["gamma"]])))
        title_parts <- if (is.null(title_parts)) gamma_part else
          bquote(paste(.(title_parts), ", ", .(gamma_part)))
      }

      if (!is.null(title_parts)) {
        graphics::title(main = title_parts)
      }
    }

    polygon_args <- utils::modifyList(
      polygon_args,
      list(
        y = c(summary_i[["hi"]], rev(summary_i[["lo"]])),
        x = c(summary_i[["alpha"]], rev(summary_i[["alpha"]]))
      )
    )
    do.call(graphics::polygon, polygon_args)

    lines_args <- utils::modifyList(
      lines_args,
      list(
        x = summary_i[["alpha"]],
        y = summary_i[["mean"]]
      )
    )
    do.call(graphics::lines, lines_args)

    if (
      plot_min &&
        optimum[["q"]] == g[["q"]] &&
        optimum[["gamma"]] == g[["gamma"]]
    ) {
      do.call(graphics::abline, abline_args)
    }
  }
  palette

  invisible()
}
