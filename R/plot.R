#' Plot coefficients
#'
#' Plot the fitted model's regression
#' coefficients along the regularization path.
#'
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap xlab ylab theme_bw theme
#' labs
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
plot.SLOPE = function(x,
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
              alpha = object$alpha,
              deviance_ratio = object$deviance_ratio,
              step = seq_along(object$alpha))

  xlab <- switch(x_variable,
                 alpha = expression(alpha),
                 deviance_ratio = "Fraction Deviance Explained",
                 step = "Step")


  d <- as.data.frame(as.table(coefs))
  d[["x"]] <- rep(x, each = p*m)

  if(m > 1) {
    ggplot(d, aes(x = x, y = Freq, col = Var1)) +
      geom_line() +
      facet_wrap(~Var2) +
      ylab(expression(hat(beta))) +
      xlab(xlab) +
      labs(color='Variable name') +
      theme_bw() +
      theme(...)
  }else {
    ggplot(d, aes(x = x, y = Freq, col = Var1)) +
      geom_line() +
      ylab(expression(hat(beta))) +
      xlab(xlab) +
      labs(color='Variable name') +
      theme_bw() +
      theme(...)
  }
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
plot.TrainedSLOPE <-
  function(x,
           measure = c("auto", "mse", "mae", "deviance", "auc", "misclass"),
           plot_min = TRUE,
           ci_alpha = 0.2,
           ci_border = FALSE,
           ci_col = lattice::trellis.par.get("superpose.line")$col,
           ...) {

  object <- x
  family <- object$model$family

  measure <- match.arg(measure)

  if (measure == "auto") {
    measure <- switch(
      family,
      gaussian = "mse",
      binomial = "deviance",
      multinomial = "deviance",
      poisson = "mse"
    )

    ind <- match(measure, object$measure$measure)

    if (is.na(ind))
      measure <- object$measure$measure[1]
  }

  ind <- match(measure, object$measure$measure)

  if (is.na(ind))
    stop("measure ", measure, " was not used or not available when",
         "fitting the model")

  if (length(measure) > 1)
    stop("you are only allowed to plot one measure at a time")

  measure_label <- object$measure$label[ind]

  summary <- object$summary[object$summary$measure == measure, ]
  optimum <- object$optima[ind, , drop = FALSE]
  model <- object$model

  alpha <- unique(summary$alpha)
  q <- unique(summary$q)

  summary$q <- as.factor(summary$q)

  # get indices of best fit
  best_ind <- match(optimum$alpha, summary$alpha)

  if (length(q) > 1) {
    x <- quote(mean ~ alpha | q)

    strip <- lattice::strip.custom(
      var.name = "q",
      sep = expression(" = "),
      strip.names = TRUE
    )
    best_outer_ind <- match(optimum$q, unique(summary$q))
  } else {
    x <- quote(mean ~ alpha)
    strip <- lattice::strip.default
    best_outer_ind <- 1
  }

  xlab <- expression(log[e](alpha))

  args <- list(
    x = x,
    data = summary,
    type = "l",
    scales = list(x = list(log = "e", relation = "free")),
    xlab = xlab,
    ylab = measure_label,
    grid = FALSE,
    lower = summary$lo,
    upper = summary$hi,
    plot_min = plot_min,

    prepanel = function(x,
                        y,
                        lower,
                        upper,
                        subscripts,
                        groups = NULL,
                        ...) {
      if (any(!is.na(x)) && any(!is.na(y))) {
        ord <- order(as.numeric(x))
        if (!is.null(groups)) {
          gg <- groups[subscripts]
          dx <- unlist(lapply(split(as.numeric(x)[ord], gg[ord]), diff))
          dy <- unlist(lapply(split(as.numeric(y)[ord], gg[ord]), diff))
        } else {
          dx <- diff(as.numeric(x[ord]))
          dy <- diff(as.numeric(y[ord]))
        }
        list(xlim = range(x, finite = TRUE),
             ylim = range(c(lower, upper), finite = TRUE),
             dx = dx,
             dy = dy,
             xat = if (is.factor(x)) sort(unique(as.numeric(x))) else NULL,
             yat = if (is.factor(y)) sort(unique(as.numeric(y))) else NULL)
      } else {
        list(xlim = rep(NA, 2),
             ylim = rep(NA, 2),
             dx = NA,
             dy = NA)
      }
    },

    xscale.components = function(lim, ...) {
      x <- lattice::xscale.components.default(lim, ...)
      x$bottom$labels$labels <- parse(text = x$bottom$labels$labels)
      x
    },

    strip = strip,

    panel = function(x,
                     y,
                     subscripts,
                     lower,
                     upper,
                     grid,
                     plot_min,
                     plot_1se,
                     ...) {
      if (isTRUE(grid))
        lattice::panel.grid(h = -1, v = -1)

      lattice::panel.polygon(
        c(x, rev(x)),
        c(upper[subscripts],
          rev(lower[subscripts])),
        col = ci_col,
        alpha = ci_alpha,
        border = ci_border
      )

      if (lattice::packet.number() == best_outer_ind) {
        if (plot_min)
          lattice::panel.refline(v = x[best_ind],
                                 col = 1,
                                 lty = 2)
      }

      lattice::panel.xyplot(x, y, ...)
    }
  )

  args <- utils::modifyList(args, list(...))

  do.call(lattice::xyplot, args)
}
