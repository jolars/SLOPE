#' Plot coefficients
#'
#' Plot the fitted model's regression
#' coefficients along the regularization path.
#'
#' @param x an object of class `"SLOPE"`
#' @param ... parameters that will be used to modify the call to
#'   [lattice::xyplot()]
#' @param intercept whether to plot the intercept
#'
#' @seealso [lattice::xyplot()], [SLOPE()], [plotDiagnostics()]
#'
#' @return An object of class `"trellis"`, which will be plotted on the
#'   current device unless stored in a variable.
#' @export
#'
#' @examples
#' fit <- SLOPE(heart$x, heart$y)
#' plot(fit)
plot.SLOPE = function(x, intercept = FALSE, ...) {
  object <- x

  coefs <- object$coefficients

  intercept_in_model <- "(Intercept)" %in% rownames(coefs)
  include_intercept <- intercept && intercept_in_model

  nz <- which(apply(object$nonzeros, 1, any))

  if (include_intercept) {
    ind <- c(1, nz + 1)
  } else {
    ind <- nz + intercept_in_model
  }

  coefs <- coefs[ind, , , drop = FALSE]

  if (is.null(coefs))
    stop("nothing to plot since model is yet to be fit")

  p <- NROW(coefs) # number of features
  m <- NCOL(coefs) # number of responses

  args <- list()

  x <- object$sigma
  xlab <- expression(sigma)

  n_x <- length(x)
  d <- as.data.frame(as.table(coefs))
  d$x <- rep(x, each = p*m)

  n_col <- length(lattice::trellis.par.get("superpose.line")$col)
  lty <- rep(1:4, each = n_col)
  superpose.line <- list(lty = lty)

  args <- list(
    x = if (m > 1)
      quote(Freq ~ x | Var2)
    else
      quote(Freq ~ x),
    type = if (n_x == 1) "p" else "l",
    groups = quote(Var1),
    data = quote(d),
    ylab = expression(hat(beta)),
    xlab = xlab,
    par.settings = list(superpose.line = superpose.line),
    # scales = list(x = list(log = "e")),
    # xscale.components = function(lim, ...) {
    #   x <- lattice::xscale.components.default(lim, ...)
    #   x$bottom$labels$labels <- parse(text = x$bottom$labels$labels)
    #   x
    # },
    auto.key = if (p <= 20)
      list(space = "right", lines = TRUE, points = FALSE)
    else FALSE,
    abline = within(lattice::trellis.par.get("reference.line"), {h = 0})
  )

  # switch(match.arg(xvar),
  #        norm = {
  #          plot_args$xlab <-
  #            expression(group("|", group("|", hat(beta), "|"), "|")[1])
  #          plot_data$xval <- if (is.list(beta))
  #            rowSums(vapply(beta,
  #                           function(x) colSums(abs(as.matrix(x))),
  #                           double(ncol(beta[[1]]))))
  #          else
  #            colSums(abs(as.matrix(beta)))
  #        },
  #        lambda = {
  #          plot_args$xlab <- expression(lambda)
  #          plot_args$scales <- list(x = list(log = "e"))
  #          plot_data$xval <- x$lambda
  #
  #          # Prettier x scale
  #          plot_args$xscale.components <- function(lim, ...) {
  #            x <- lattice::xscale.components.default(lim, ...)
  #            x$bottom$labels$labels <- parse(text = x$bottom$labels$labels)
  #            x
  #          }
  #        },
  #        dev = {
  #          plot_args$xlab <- "Fraction of deviance explained"
  #          plot_data$xval <- x$dev.ratio
  #        })

  # Let the user modify the plot parameters
  do.call(lattice::xyplot, utils::modifyList(args, list(...)))
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

  sigma <- unique(summary$sigma)
  q <- unique(summary$q)

  summary$q <- as.factor(summary$q)

  # get indices of best fit
  best_ind <- match(optimum$sigma, summary$sigma)

  if (length(q) > 1) {
    x <- quote(mean ~ sigma | q)

    strip <- lattice::strip.custom(
      var.name = "q",
      sep = expression(" = "),
      strip.names = TRUE
    )
    best_outer_ind <- match(optimum$q, unique(summary$q))
  } else {
    x <- quote(mean ~ sigma)
    strip <- lattice::strip.default
    best_outer_ind <- 1
  }

  xlab <- expression(log[e](sigma))

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
