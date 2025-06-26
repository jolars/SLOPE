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
#' @param add_labels whether to add labels (numbers) on the right side
#'   of the plot for each coefficient
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
  add_labels = FALSE,
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
    col = grDevices::palette.colors(10, "Tableau")
  )

  plot_args <- utils::modifyList(
    default_plot_args,
    list(...)
  )

  if (m == 1) {
    coefs_k <- t(coefs)
    xlim <- if (x_variable == "alpha") rev(range(x)) else range(x)
    plot_args_k <- utils::modifyList(plot_args, list(x = x, y = coefs_k))
    do.call(graphics::matplot, plot_args_k)

    if (add_labels) {
      addCoefLabels(coefs_k, x)
    }
  } else {
    for (k in seq_len(m)) {
      coefs_k <- t(coefs[[k]])
      plot_args_k <- utils::modifyList(
        plot_args,
        list(
          x = x,
          y = coefs_k,
          main = paste0("Response: ", k)
        )
      )
      do.call(graphics::matplot, plot_args_k)

      if (add_labels) {
        addCoefLabels(coefs_k, x)
      }
    }
  }

  invisible()
}

addCoefLabels <- function(coef_matrix, x) {
  n_coefs <- ncol(coef_matrix)
  last_x <- x[length(x)]

  for (i in seq_len(n_coefs)) {
    last_y <- coef_matrix[nrow(coef_matrix), i]
    graphics::text(last_x, last_y, i, pos = 4, cex = 0.8)
  }
}

#' Plot results from cross-validation
#'
#' @param x an object of class `'TrainedSLOPE'`, typically from a call
#'   to [cvSLOPE()]
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
#' @param index an optional index, to plot only one (the index-th) set
#'   of the parameter combinations.
#' @param ... ignored
#'
#' @seealso [cvSLOPE()]
#' @family model-tuning
#'
#' @return A plot for every value of `q` is produced on the current device.
#'
#' @export
#'
#' @examples
#' # Cross-validation for a SLOPE binomial model
#' set.seed(123)
#' tune <- cvSLOPE(
#'   subset(mtcars, select = c("mpg", "drat", "wt")),
#'   mtcars$hp,
#'   q = c(0.1, 0.2),
#'   n_folds = 10
#' )
#' plot(tune, ci_col = "salmon", index = 1)
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
  index = NULL,
  measure,
  ...
) {
  object <- x

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

  multi_param <- nrow(grid) > 1

  if (!is.null(index)) {
    stopifnot(
      is.numeric(index),
      index >= 1,
      index <= nrow(grid)
    )
    grid <- grid[index, , drop = FALSE]
  }

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]

    ind <- which(
      summary[["q"]] == g[["q"]] &
        summary[["gamma"]] == g[["gamma"]]
    )
    summary_i <- summary[ind, ]

    # Plot frame
    x <- summary_i[["alpha"]]
    plot_args <- utils::modifyList(
      plot_args,
      list(
        x = x,
        xlim = rev(range(x)),
        y = summary_i[["mean"]]
      )
    )

    do.call(plot, plot_args)

    if (multi_param) {
      title_parts <- NULL
      if (length(q) > 1) {
        q_part <- bquote(paste("q = ", .(g[["q"]])))
        title_parts <- if (is.null(title_parts)) {
          q_part
        } else {
          bquote(paste(.(title_parts), ", ", .(q_part)))
        }
      }
      if (length(gamma) > 1) {
        gamma_part <- bquote(paste(gamma, " = ", .(g[["gamma"]])))
        title_parts <- if (is.null(title_parts)) {
          gamma_part
        } else {
          bquote(paste(.(title_parts), ", ", .(gamma_part)))
        }
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

  invisible()
}


#' Plot cluster structure
#'
#' Note that this function requires the `patterns` argument to be set to
#' `TRUE` in the call to [SLOPE()]. Calling this function on a
#' `SLOPE` object without patterns will result in an error.
#'
#' @param x an object of class `'SLOPE'`
#' @param plot_signs logical, indicating whether to plot signs of estimated
#' coefficients on the plot
#' @param color_clusters logical, indicating whether the clusters should have
#' different colors
#' @param include_zeroes logical, indicating whether zero variables should be
#' plotted. Default to TRUE
#' @param show_alpha logical, indicating whether labels with alpha values or
#' steps in the path should be plotted.
#' @param alpha_steps a vector of integer alpha steps to plot. If `NULL`,
#' all the steps are plotted.
#' @param palette a character string specifying the color palette to use for
#' the clusters. This is passed to [grDevices::hcl.colors()].
#' @param ... additional arguments passed to [graphics::image()].
#'
#' @seealso [SLOPE()], [graphics::image()], [graphics::text()].
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
#' plotClusters(fit)
#' plotClusters(fit, alpha_steps = 1:10)
plotClusters <- function(
  x,
  plot_signs = FALSE,
  color_clusters = TRUE,
  include_zeroes = TRUE,
  show_alpha = FALSE,
  alpha_steps = NULL,
  palette = "viridis",
  ...
) {
  object <- x

  pat <- object$patterns

  if (is.null(pat)) {
    stop(
      "No patterns found in the SLOPE object.",
      "Please run SLOPE with `patterns = TRUE`."
    )
  }

  pat[[1]] <- as.matrix(numeric(length(object$lambda)), ncol = 1)

  if (is.null(alpha_steps)) {
    alpha_steps <- seq_along(pat)
  } else {
    alpha_steps <- unique(alpha_steps)
    stopifnot(
      is.numeric(alpha_steps),
      all(alpha_steps %% 1 == 0),
      all(alpha_steps >= 1),
      all(alpha_steps <= length(pat))
    )
  }

  mat <- sapply(pat[alpha_steps], function(m) {
    rowSums(t(t(as.matrix(m)) * seq_len(ncol(as.matrix(m)))))
  })

  rownames(mat) <- seq_len(nrow(mat))

  if (!include_zeroes) {
    mat <- mat[rowSums(mat) != 0, ]
  }

  abs_mat <- abs(mat)

  abs_vals <- sort(unique(as.vector(abs_mat)))

  if (color_clusters) {
    n_clusters <- length(abs_vals)

    my_colors <- c(
      "white",
      grDevices::hcl.colors(
        n_clusters - 1,
        palette = palette,
        alpha = NULL,
        rev = FALSE,
        fixup = TRUE
      )
    )
  } else {
    my_colors <- c("white", rep("grey", length(abs_vals) - 1))
  }

  if (show_alpha) {
    step <- round(object$alpha, 3)
    xlabel <- expression(alpha)
  } else {
    step <- seq_len(ncol(mat))
    xlabel <- "Step"
  }

  breaks <- c(-1, abs_vals)

  default_image_args <- list(
    col = my_colors,
    breaks = breaks,
    axes = FALSE,
    xlab = xlabel,
    ylab = "Variable"
  )

  image_args <- utils::modifyList(default_image_args, list(...))

  do.call(graphics::image, c(list(t(abs_mat)), image_args))

  x_coords <- seq(0, 1, length.out = ncol(mat))
  x_coords <- x_coords + mean(x_coords[1:2])

  y_coords <- seq(0, 1, length.out = nrow(mat))
  y_coords <- y_coords + mean(y_coords[1:2])

  if (plot_signs) {
    for (i in seq_len(nrow(mat))) {
      for (j in seq_len(ncol(mat))) {
        val <- mat[rev(seq_len(nrow(mat))), ][i, j]
        sign_char <- ifelse(val > 0, "+", ifelse(val < 0, "-", ""))

        x <- (j - 1) / (ncol(mat) - 1)
        y <- 1 - (i - 1) / (nrow(mat) - 1)

        graphics::text(x, y, labels = sign_char)
      }
    }
  }

  graphics::abline(v = x_coords, col = "white")
  graphics::abline(h = y_coords, col = "white")

  graphics::box()

  graphics::axis(
    1,
    at = seq(0, 1, length.out = ncol(mat)),
    labels = step,
    las = 2
  )
  graphics::axis(
    2,
    at = seq(0, 1, length.out = nrow(mat)),
    labels = rownames(mat),
    las = 1
  )

  invisible()
}
