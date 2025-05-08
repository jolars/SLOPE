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
  xlim <- if (x_variable == "alpha") rev(range(x)) else range(x)
  log_var <- if (x_variable == "alpha") "x" else ""

  if (m == 1) {
    xlim <- if (x_variable == "alpha") rev(range(x)) else range(x)
    plot_coefs(x, coefs, xlab, xlim, log_var, ...)
  } else {
    for (k in seq_len(m)) {
      plot_coefs(
        x,
        coefs[[k]],
        xlab,
        xlim,
        log_var,
        main = paste0("Response: ", k),
        ...
      )
    }
  }

  invisible()
}


plot_coefs <- function(x, coefs, xlab, xlim, log_var = "x", ...) {
  graphics::matplot(
    x,
    t(coefs),
    type = "l",
    lty = 1,
    xlab = xlab,
    xlim = xlim,
    ylab = expression(hat(beta)),
    log = log_var,
    ...
  )
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
  ci_border = NA,
  ci_col = "salmon",
  plot_args = list(),
  polygon_args = list(),
  lines_args = list(),
  abline_args = list(),
  ...
) {
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

  for (i in seq_along(q)) {
    summary_i <- summary[summary[["q"]] == q[i], ]

    # Plot frame
    plot_args <- utils::modifyList(
      plot_args,
      list(
        x = summary_i[["alpha"]],
        y = summary_i[["mean"]]
      )
    )
    do.call(plot, plot_args)

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

    if (plot_min && optimum[["q"]] == q[i]) {
      do.call(graphics::abline, abline_args)
    }
  }

  invisible()
}



#' Plot cluster structure
#'
#' @param x an object of class `'SLOPE'`
#' @param plot_signs logical, indicating whether to plot signs of estimated
#' coefficients on the plot
#' @param color_clusters logical, indicating whether the clusters should have
#' different colors
#' @param include_zeroes logical, indicating whether zero variables should be
#' plotted. Default to TRUE
#' @param show_alpha logical, indicatiung whether labels with alpha values or
#' steps in the path should be plotted. Default FALSE.
#' @param alpha_steps a vector of integer alpha steps to plot. If NULL, all the steps
#' are plotted. Default to NULL.
#'
#' @seealso [SLOPE()]
#'
#' @return Invisibly returns NULL. The function is called for its
#'   side effect of producing a plot.
#'
#' @export
#'
#' @examples
#' set.seed(10)
#' X <- matrix(rnorm(10000), ncol = 10)
#' colnames(X) <- paste0("X", 1:10)
#' beta <- c(rep(10, 3), rep(-20, 2), rep(20, 2), rep(0, 3))
#' Y <- X %*% beta + rnorm(1000)
#' fit <- SLOPE(X, Y, patterns = TRUE)
#'
#' plot_clusters(fit)
#' plot_clusters(fit, alpha_steps = 1:10)

plot_clusters <- function(
    x, plot_signs = FALSE,
    color_clusters = TRUE,
    include_zeroes = TRUE,
    show_alpha = FALSE,
    alpha_steps = NULL
) {
  object <- x

  pat <- object$patterns
  pat[[1]] <- as.matrix(numeric(length(object$lambda)), ncol = 1)

  if (is.null(alpha_steps)) {
    alpha_steps <- 1:length(pat)
  } else {
    alpha_steps <- unique(alpha_steps)
    stopifnot(is.numeric(alpha_steps),
              all(alpha_steps %% 1 == 0),
              all(alpha_steps >= 1),
              all(alpha_steps <= length(pat)))
  }

  mat <- sapply(pat[alpha_steps], function(m) {
    rowSums(t(t(as.matrix(m)) * 1:ncol(as.matrix(m))))
  })

  rownames(mat) <- 1:nrow(mat)

  if (!include_zeroes) mat <- mat[rowSums(mat) != 0, ]

  abs_mat <- abs(mat)

  abs_vals <- sort(unique(as.vector(abs_mat)))

  if (color_clusters){
    my_colors <- c("white", grDevices::rainbow(length(abs_vals) - 1, alpha = 0.7))
  } else {
    my_colors <- c("white", rep("grey", length(abs_vals) - 1))
  }

  if (show_alpha) {
    step <- round(object$alpha, 3)
    xlabel <- "alpha"
  } else {
    step <- 1:ncol(mat)
    xlabel <- "path step"
  }

  breaks <- c(-1, abs_vals)

  graphics::image(t(abs_mat),
                  col = my_colors,
                  breaks = breaks,
                  axes = FALSE,
                  xlab = xlabel, ylab = "variable")

  graphics::box(col = "black", lwd = 1.5)
  graphics::axis(1,
                 at = seq(0, 1, length.out = ncol(mat)),
                 labels = step,
                 las = 2)
  graphics::axis(2,
                 at = seq(0, 1, length.out = nrow(mat)),
                 labels = rownames(mat),
                 las = 1)

  if (plot_signs) {
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        val <- mat[nrow(mat):1, ][i, j]
        sign_char <- ifelse(val > 0, "+", ifelse(val < 0, "-", ""))

        x <- (j - 1) / (ncol(mat) - 1)
        y <- 1 - (i - 1) / (nrow(mat) - 1)

        graphics::text(x, y, labels = sign_char)
      }
    }
  }

  x_coords <- seq(0, 1, length.out = ncol(mat))
  x_coords <- x_coords + mean(x_coords[1:2])

  y_coords <- seq(0, 1, length.out = nrow(mat))
  y_coords <- y_coords + mean(y_coords[1:2])

  graphics::abline(v = x_coords, col = grDevices::adjustcolor("black", alpha = 0.5), lwd = 0.5)
  graphics::abline(h = y_coords, col = grDevices::adjustcolor("black", alpha = 0.5), lwd = 0.5)

  invisible()
}
