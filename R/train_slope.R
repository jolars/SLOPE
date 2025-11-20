#' Train a SLOPE Model
#'
#' This function trains a model fit by [SLOPE()] by tuning its parameters
#' through cross-validation.
#'
#' Note that by default this method matches all of the available metrics
#' for the given model family against those provided in the argument
#' `measure`. Collecting these measures is not particularly demanding
#' computationally so it is almost always best to leave this argument
#' as it is and then choose which argument to focus on in the call
#' to [plot.TrainedSLOPE()].
#'
#' @inheritParams SLOPE
#' @param number number of folds (cross-validation)
#' @param repeats number of repeats for each fold (for repeated *k*-fold cross
#'   validation)
#' @param measure measure to try to optimize; note that you may
#'   supply *multiple* values here and that, by default,
#'   all the possible measures for the given model will be used.
#' @param ... other arguments to pass on to [SLOPE()]
#'
#' @return An object of class `"TrainedSLOPE"`, with the following slots:
#' \item{summary}{a summary of the results with means, standard errors,
#'                and 0.95 confidence levels}
#' \item{data}{the raw data from the model training}
#' \item{optima}{a `data.frame` of the best (mean)
#'   values for the different metrics and their corresponding parameter values}
#' \item{measure}{a `data.frame` listing the used metrics and their labels}
#' \item{model}{the model fit to the entire data set}
#' \item{call}{the call}
#'
#' @export
#' @family model-tuning
#'
#' @examples
#' # 8-fold cross-validation repeated 5 times
#' tune <- trainSLOPE(subset(mtcars, select = c("mpg", "drat", "wt")),
#'   mtcars$hp,
#'   q = c(0.1, 0.2),
#'   number = 8,
#'   repeats = 5,
#'   measure = "mse"
#' )
trainSLOPE <- function(
  x,
  y,
  q = 0.2,
  number = 10,
  repeats = 1,
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

  if ("missclass" %in% measure) {
    warning("measure 'missclass' is deprecated; please use 'misclass' instead.")

    measure <- measure[!(measure %in% "missclass")]

    if (!("misclass" %in% measure)) {
      measure <- c(measure, "misclass")
    }
  }

  n <- NROW(x)

  measure <- match.arg(measure, several.ok = TRUE)

  y <- as.matrix(y)

  stopifnot(
    NROW(x) > number,
    number > 1,
    repeats >= 1
  )

  # get initial penalty sequence
  fit <- SLOPE(x, y, ...)

  # match measure against accepted measure for the given family
  family <- fit$family

  ok <- switch(
    family,
    gaussian = c("mse", "mae"),
    binomial = c("mse", "mae", "deviance", "misclass", "auc"),
    poisson = c("mse", "mae"),
    multinomial = c("mse", "mae", "deviance", "misclass")
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

  alpha <- fit$alpha
  path_length <- length(alpha)
  n_q <- length(q)
  n_measure <- length(measure)

  fold_size <- ceiling(n / number)

  # list of repeated folds
  fold_id <- rep(
    list(matrix(
      c(sample(n), rep(0, number * fold_size - n)),
      fold_size,
      byrow = TRUE
    )),
    repeats
  )

  grid <- expand.grid(
    q = q,
    fold = seq_len(number),
    repetition = seq_len(repeats)
  )

  # Initialize list to store results
  r <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    id <- grid[["fold"]][i]
    repetition <- grid[["repetition"]][i]
    q <- grid[["q"]][i]

    test_ind <- fold_id[[repetition]][, id]

    x_train <- x[-test_ind, , drop = FALSE]
    y_train <- y[-test_ind, , drop = FALSE]
    x_test <- x[test_ind, , drop = FALSE]
    y_test <- y[test_ind, , drop = FALSE]

    # arguments for SLOPE
    args <- utils::modifyList(
      list(
        x = x_train,
        y = y_train,
        q = q,
        alpha = alpha
      ),
      list(...)
    )

    # fitting model
    fit_id <- do.call(SLOPE::SLOPE, args)

    s <- lapply(measure, function(m) {
      SLOPE::score(fit_id, x_test, y_test, measure = m)
    })

    r[[i]] <- unlist(s)
  }

  tmp <- array(unlist(r), c(path_length * n_q, n_measure, number * repeats))
  d <- matrix(tmp, c(path_length * n_q * n_measure, number * repeats))

  means <- rowMeans(d)
  se <- apply(d, 1, stats::sd) / sqrt(repeats * number)
  ci <- stats::qt(0.975, number * repeats - 1) * se
  lo <- means - ci
  hi <- means + ci

  summary <- data.frame(
    q = rep(q, each = path_length * n_measure),
    alpha = rep(alpha, n_measure * n_q),
    measure = rep(measure, each = path_length, times = n_q),
    mean = means,
    se = se,
    lo = lo,
    hi = hi
  )

  optima <- do.call(
    rbind,
    by(
      summary,
      as.factor(summary$measure),
      function(x) {
        if (x$measure[1] == "auc") {
          x[which.max(x$mean), ]
        } else {
          x[which.min(x$mean), ]
        }
      }
    )
  )

  labels <- vapply(
    measure,
    function(m) {
      switch(
        m,
        deviance = {
          if (inherits(fit, "GaussianSLOPE")) {
            "Mean-Squared Error"
          } else if (inherits(fit, "BinomialSLOPE")) {
            "Binomial Deviance"
          } else if (inherits(fit, "PoissonSLOPE")) {
            "Mean-Squared Error"
          } else if (inherits(fit, "MultinomialSLOPE")) {
            "Multinomial Deviance"
          }
        },
        mse = "Mean Squared Error",
        mae = "Mean Absolute Error",
        accuracy = "Accuracy",
        auc = "AUC",
        misclass = "Misclassification Rate"
      )
    },
    FUN.VALUE = character(1)
  )

  rownames(summary) <- NULL
  rownames(optima) <- NULL

  structure(
    list(
      summary = summary,
      data = d,
      optima = optima,
      measure = data.frame(
        measure = measure,
        label = labels,
        row.names = NULL,
        stringsAsFactors = FALSE
      ),
      model = fit,
      call = ocall
    ),
    class = "TrainedSLOPE"
  )
}
