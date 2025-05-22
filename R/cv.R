#' Tune SLOPE with cross-validation
#'
#' This function trains a model fit by [SLOPE()] by tuning its parameters
#' through cross-validation.
#'
#' @inheritParams SLOPE
#' @param q a vector of quantiles for the `q` parameter in SLOPE
#' @param gamma relaxation parameter for SLOPE. Default is `0.0`, which
#'   implies to relaxation of the penalty.
#' @param n_folds number of folds (cross-validation)
#' @param n_repeats number of folds (cross-validation)
#' @param measure DEPRECATED
#' @param ... other arguments to pass on to [SLOPE()]
#'
#' @return An object of class `"TrainedSLOPE"`, with the following slots:
#' \item{summary}{a summary of the results with means, standard errors,
#'                and 0.95 confidence levels}
#' \item{data}{the raw data from the model training}
#' \item{optima}{a `data.frame` of the best (mean)
#'   values for the different metrics and their corresponding parameter values}
#' \item{measure}{a `data.frame` listing the used metric and its label}
#' \item{call}{the call}
#'
#' @export
#'
#' @seealso [plot.TrainedSLOPE()]
#' @family model-tuning
#'
#' @examples
#' # 8-fold cross-validation
#' tune <- cvSLOPE(
#'   subset(mtcars, select = c("mpg", "drat", "wt")),
#'   mtcars$hp,
#'   q = c(0.1, 0.2),
#'   n_folds = 8,
#'   n_repeats = 2,
#'   measure = "mse"
#' )
cvSLOPE <- function(
  x,
  y,
  q = 0.2,
  gamma = 0.0,
  n_folds = 10,
  n_repeats = 1,
  measure = c(
    "mse",
    "mae",
    "deviance",
    "misclass",
    "auc"
  ),
  ...
) {
  ocall <- match.call()

  n <- NROW(x)

  measure <- match.arg(measure)

  stopifnot(
    NROW(x) > n_folds,
    n_folds > 1
  )

  slope_args <- utils::modifyList(list(...), list(x = x, y = y))
  model_args <- do.call(processSlopeArgs, slope_args)

  x <- model_args$x
  y <- model_args$y

  family <- model_args[["family"]]

  ok <- switch(
    family,
    gaussian = c("deviance", "mse", "mae"),
    binomial = c("deviance", "mse", "mae", "misclass", "auc"),
    poisson = c("deviance", "mse", "mae"),
    multinomial = c("deviance", "mse", "mae", "misclass", "auc")
  )

  measure <- measure[measure %in% ok]

  if (length(measure) == 0) {
    stop(paste0(
      "For the given family: ",
      family,
      ", measure needs to be one of: ",
      paste0(ok, collapse = ", ")
    ))
  }

  folds <- createFolds(n, n_folds, n_repeats)

  cvCpp <- if (inherits(x, "sparseMatrix")) cvSparseCpp else cvDenseCpp

  cv_args <- list(
    q = q,
    gamma = gamma,
    metric = measure,
    predefined_folds = folds
  )

  res <- cvCpp(x, y, cv_args, model_args)

  summary_list <- lapply(res[["results"]], function(results_i) {
    means <- results_i[["mean_scores"]]
    se <- results_i[["std_errors"]]
    ci <- stats::qt(0.975, NROW(results_i[["score"]]) - 1) * se

    data.frame(
      q = results_i[["params"]][["q"]],
      gamma = results_i[["params"]][["gamma"]],
      alpha = results_i[["alphas"]],
      measure = measure,
      mean = means,
      se = se,
      lo = means - ci,
      hi = means + ci
    )
  })
  summary <- do.call(rbind, summary_list)

  comparator <- if (measure %in% c("deviance", "mse", "mae")) {
    which.min
  } else {
    which.max
  }

  optima <- summary[comparator(summary$mean), , drop = FALSE]

  labels <- vapply(
    measure,
    function(m) {
      switch(
        m,
        deviance = switch(
          family,
          gaussian = "Mean-Squared Error",
          binomial = "Deviance",
          poisson = "Deviance",
          multinomial = "Deviance"
        ),
        mse = "Mean Squared Error",
        mae = "Mean Absolute Error",
        accuracy = "Accuracy",
        auc = "AUC",
        misclass = "Misclassification Rate"
      )
    },
    FUN.VALUE = character(1)
  )

  structure(
    list(
      summary = summary,
      data = res$results,
      optima = optima,
      measure = data.frame(
        measure = measure,
        label = labels,
        row.names = NULL
      ),
      call = ocall
    ),
    class = "TrainedSLOPE"
  )
}

#' Create cross-validation folds
#'
#' Internal function that creates fold assignments for k-fold cross-validation,
#' with support for repeated cross-validation.
#'
#' @param n Integer. Number of observations to split into folds.
#' @param n_folds Integer. Number of folds to create.
#' @param n_repeats Integer. Number of times to repeat the cross-validation
#'   process with different fold assignments. Default is 1.
#'
#' @return A list of length `n_repeats`. Each element contains a list of
#'   `n_folds` integer vectors representing the indices (0-based)
#'   of observations in each fold.
#'
#' @keywords internal
createFolds <- function(n, n_folds, n_repeats = 1) {
  repeated_folds <- vector("list", n_repeats)

  for (i in seq_len(n_repeats)) {
    fold_assignments <- rep(seq_len(n_folds), length.out = n)
    fold_assignments <- sample(fold_assignments)
    folds <- split(seq_len(n) - 1, fold_assignments)
    names(folds) <- NULL

    repeated_folds[[i]] <- folds
  }

  repeated_folds
}
