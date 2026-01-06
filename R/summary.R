#' Summarize SLOPE Model
#'
#' Produces a summary of a fitted SLOPE model, including information about
#' the regularization path, model family, and fitted values.
#'
#' @param object an object of class `"SLOPE"`, typically from a call to
#'   [SLOPE()]
#' @param ... other arguments (currently ignored)
#'
#' @return An object of class `"summary_SLOPE"` with the following components:
#' \item{call}{the call that produced the model}
#' \item{family}{the model family}
#' \item{n_obs}{number of observations}
#' \item{n_predictors}{number of predictors}
#' \item{has_intercept}{whether an intercept was fit}
#' \item{path_length}{number of steps in the regularization path}
#' \item{alpha_range}{range of alpha values in the path}
#' \item{deviance_ratio_range}{range of deviance ratios in the path}
#' \item{null_deviance}{null deviance}
#' \item{path_summary}{data frame summarizing the regularization path}
#'
#' @seealso [SLOPE()], [print.summary_SLOPE()]
#' @family SLOPE-methods
#'
#' @examples
#' fit <- SLOPE(heart$x, heart$y)
#' summary(fit)
#'
#' # Multinomial example
#' fit_multi <- SLOPE(wine$x, wine$y, family = "multinomial")
#' summary(fit_multi)
#'
#' @method summary SLOPE
#' @export
summary.SLOPE <- function(object, ...) {
  n_obs <- object$n_observations
  n_predictors <- object$n_predictors

  path_length <- length(object$alpha)
  alpha_range <- range(object$alpha)
  deviance_ratio_range <- range(object$deviance_ratio)

  n_nonzero <- vapply(object$nonzeros, sum, FUN.VALUE = double(1))

  path_summary <- data.frame(
    alpha = object$alpha,
    deviance_ratio = object$deviance_ratio,
    n_nonzero = n_nonzero
  )

  structure(
    list(
      call = object$call,
      family = object$family,
      n_obs = n_obs,
      n_predictors = n_predictors,
      has_intercept = object$has_intercept,
      path_length = path_length,
      alpha_range = alpha_range,
      deviance_ratio_range = deviance_ratio_range,
      null_deviance = object$null_deviance,
      path_summary = path_summary
    ),
    class = "summary_SLOPE"
  )
}

#' Print Summary of SLOPE Model
#'
#' @param x an object of class `'summary_SLOPE'`
#' @param digits number of significant digits to print
#' @param ... other arguments passed to [print()]
#'
#' @return Invisibly returns the input object
#'
#' @seealso [summary.SLOPE()]
#' @method print summary_SLOPE
#' @export
print.summary_SLOPE <- function(x, digits = 3, ...) {
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n")

  cat("Family:", x$family, "\n")
  if (!is.na(x$n_obs)) {
    cat("Observations:", x$n_obs, "\n")
  }
  cat("Predictors:", x$n_predictors, "\n")
  cat("Intercept:", if (x$has_intercept) "Yes" else "No", "\n")
  cat("\nRegularization path:\n")
  cat("  Length:", x$path_length, "steps\n")
  cat(
    "  Alpha range:",
    format(x$alpha_range[1], digits = digits),
    "to",
    format(x$alpha_range[2], digits = digits),
    "\n"
  )
  cat(
    "  Deviance ratio range:",
    format(x$deviance_ratio_range[1], digits = digits),
    "to",
    format(x$deviance_ratio_range[2], digits = digits),
    "\n"
  )
  cat("  Null deviance:", format(x$null_deviance, digits = digits), "\n")

  cat("\nPath summary (first and last 5 steps):\n")
  n_rows <- nrow(x$path_summary)

  # Format the path summary for better display
  formatted_summary <- data.frame(
    alpha = format(signif(x$path_summary$alpha, digits), scientific = FALSE),
    deviance_ratio = format(signif(x$path_summary$deviance_ratio, digits)),
    n_nonzero = x$path_summary$n_nonzero
  )

  if (n_rows <= 10) {
    print(formatted_summary, row.names = FALSE, ...)
  } else {
    # Print header and first 5 rows
    print(head(formatted_summary, 5), row.names = FALSE, ...)
    cat(" . . .\n")
    # Print last 5 rows without header
    tail_rows <- tail(formatted_summary, 5)
    for (i in seq_len(nrow(tail_rows))) {
      cat(sprintf(
        "%9s %14s %9d\n",
        tail_rows$alpha[i],
        tail_rows$deviance_ratio[i],
        tail_rows$n_nonzero[i]
      ))
    }
  }

  invisible(x)
}

#' Summarize TrainedSLOPE Model
#'
#' Produces a summary of a trained SLOPE model from cross-validation,
#' including information about the optimal parameters and performance metrics.
#'
#' @param object an object of class `'TrainedSLOPE'`, typically from a call to
#'   [cvSLOPE()] or [trainSLOPE()]
#' @param ... other arguments (currently ignored)
#'
#' @return An object of class `'summary_TrainedSLOPE'` with the following
#'   components:
#' \item{call}{the call that produced the model}
#' \item{measure}{the performance measure(s) used}
#' \item{optima}{optimal parameter values and corresponding performance}
#' \item{n_folds}{number of cross-validation folds}
#' \item{n_repeats}{number of cross-validation repeats}
#' \item{n_models}{total number of models evaluated}
#'
#' @seealso [cvSLOPE()], [trainSLOPE()], [print.summary_TrainedSLOPE()]
#' @family model-tuning
#'
#' @examples
#' tune <- cvSLOPE(
#'   subset(mtcars, select = c("mpg", "drat", "wt")),
#'   mtcars$hp,
#'   q = c(0.1, 0.2),
#'   n_folds = 5
#' )
#' summary(tune)
#'
#' @method summary TrainedSLOPE
#' @export
summary.TrainedSLOPE <- function(object, ...) {
  n_models <- nrow(object$summary)

  # Try to extract fold information from the call
  call_list <- as.list(object$call)
  n_folds <- if ("n_folds" %in% names(call_list)) {
    eval(call_list$n_folds)
  } else {
    10 # default value
  }
  n_repeats <- if ("n_repeats" %in% names(call_list)) {
    eval(call_list$n_repeats)
  } else {
    1 # default value
  }

  structure(
    list(
      call = object$call,
      measure = object$measure,
      optima = object$optima,
      n_folds = n_folds,
      n_repeats = n_repeats,
      n_models = n_models
    ),
    class = "summary_TrainedSLOPE"
  )
}

#' Print Summary of TrainedSLOPE Model
#'
#' @param x an object of class `'summary_TrainedSLOPE'`
#' @param digits number of significant digits to print
#' @param ... other arguments passed to [print()]
#'
#' @return Invisibly returns the input object
#'
#' @seealso [summary.TrainedSLOPE()]
#' @method print summary_TrainedSLOPE
#' @export
print.summary_TrainedSLOPE <- function(x, digits = 3, ...) {
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n")

  cat("Cross-validation:\n")
  if (!is.na(x$n_folds)) {
    cat("  Folds:", x$n_folds, "\n")
  }
  if (!is.na(x$n_repeats)) {
    cat("  Repeats:", x$n_repeats, "\n")
  }
  cat("  Models evaluated:", x$n_models, "\n")

  cat("\nPerformance measure:", x$measure$label, "\n")

  cat("\nOptimal parameters:\n")

  # Format the optima data frame for display
  formatted_optima <- x$optima
  if (is.numeric(formatted_optima$q)) {
    formatted_optima$q <- signif(formatted_optima$q, digits)
  }
  if (is.numeric(formatted_optima$gamma)) {
    formatted_optima$gamma <- signif(formatted_optima$gamma, digits)
  }
  if (is.numeric(formatted_optima$alpha)) {
    formatted_optima$alpha <- signif(formatted_optima$alpha, digits)
  }
  if (is.numeric(formatted_optima$mean)) {
    formatted_optima$mean <- signif(formatted_optima$mean, digits)
  }
  if (is.numeric(formatted_optima$se)) {
    formatted_optima$se <- signif(formatted_optima$se, digits)
  }
  if (is.numeric(formatted_optima$lo)) {
    formatted_optima$lo <- signif(formatted_optima$lo, digits)
  }
  if (is.numeric(formatted_optima$hi)) {
    formatted_optima$hi <- signif(formatted_optima$hi, digits)
  }

  print(formatted_optima, row.names = FALSE, ...)

  invisible(x)
}
