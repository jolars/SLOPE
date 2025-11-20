#' Compute One of Several Loss Metrics on a New Data Set
#'
#' This function is a unified interface to return various types of loss for a
#' model fit with [SLOPE()].
#'
#' @param object an object of class `"SLOPE"`
#' @param x feature matrix
#' @param y response
#' @param measure type of target measure. `"mse"` returns mean squared error.
#'   `"mae"` returns mean absolute error, `"misclass"` returns
#'   misclassification rate, and `"auc"` returns area under the ROC curve.
#'
#' @return The measure along the regularization path depending on the
#'   value in `measure`.#'
#'
#' @export
#'
#' @family SLOPE-methods
#' @seealso [SLOPE()], [predict.SLOPE()]
#'
#' @examples
#' x <- subset(infert, select = c("induced", "age", "pooled.stratum"))
#' y <- infert$case
#'
#' fit <- SLOPE(x, y, family = "binomial")
#' score(fit, x, y, measure = "auc")
score <- function(object, x, y, measure) {
  UseMethod("score")
}

#' @rdname score
#' @export
score.GaussianSLOPE <- function(object, x, y, measure = c("mse", "mae")) {
  measure <- match.arg(measure)

  y <- as.vector(y)
  y_hat <- stats::predict(object, x, simplify = FALSE)

  switch(
    measure,
    mse = apply((y_hat - y)^2, 3, mean),
    mae = apply(abs(y_hat - y), 3, mean)
  )
}

#' @rdname score
#' @export
score.BinomialSLOPE <- function(
  object,
  x,
  y,
  measure = c(
    "mse",
    "mae",
    "deviance",
    "misclass",
    "auc"
  )
) {
  measure <- match.arg(measure)

  prob_min <- 1e-05
  prob_max <- 1 - prob_min

  y <- as.factor(y)
  y <- diag(2)[as.numeric(y), ]

  y_hat <- stats::predict(object, x, type = "response", simplify = FALSE)

  switch(
    measure,
    auc = apply(y_hat, 3, function(y_hat_i) auc(y, y_hat_i)),
    mse = apply((y_hat + y[, 1] - 1)^2 + (y_hat - y[, 2])^2, 3, mean),
    mae = apply(abs(y_hat + y[, 1] - 1) + abs(y_hat - y[, 2]), 3, mean),
    deviance = {
      y_hat <- pmin(pmax(y_hat, prob_min), prob_max)
      lp <- y[, 1] * log(1 - y_hat) + y[, 2] * log(y_hat)
      ly <- log(y)
      ly[y == 0] <- 0
      ly <- drop((y * ly) %*% c(1, 1))
      apply(2 * (ly - lp), 3, mean)
    },
    misclass = apply(y[, 1] * (y_hat > 0.5) + y[, 2] * (y_hat <= 0.5), 3, mean)
  )
}

#' @rdname score
#' @export
score.MultinomialSLOPE <- function(
  object,
  x,
  y,
  measure = c(
    "mse",
    "mae",
    "deviance",
    "misclass"
  )
) {
  measure <- match.arg(measure)
  prob_min <- 1e-05
  prob_max <- 1 - prob_min

  y_observed <- y
  y_hat <- stats::predict(object, x, type = "response", simplify = FALSE)
  y_hat_class <- stats::predict(object, x, type = "class", simplify = FALSE)

  y <- as.factor(y)
  n_levels <- length(unique(y))

  y <- diag(n_levels)[as.numeric(y), ]

  switch(
    measure,
    mse = apply(y_hat, 3, function(x) mean((x - y)^2)),
    mae = apply(y_hat, 3, function(x) mean(abs(x - y))),
    deviance = {
      y_hat <- pmin(pmax(y_hat, prob_min), prob_max)

      apply(y_hat, 3, function(p_hat) {
        -2 * sum(y * log(p_hat))
      })
    },
    misclass = apply(y_hat, 3, function(x) {
      n_obs <- nrow(x)
      x_class <- array(0, dim = dim(x))
      x_class[cbind(1:n_obs, apply(x, 1, which.max))] <- 1

      1 - sum(apply(x_class == y, 1, all)) / n_obs
    })
  )
}

#' @rdname score
#' @export
score.PoissonSLOPE <- function(object, x, y, measure = c("mse", "mae")) {
  measure <- match.arg(measure)

  y <- as.vector(y)
  y_hat <- stats::predict(object, x, simplify = FALSE)

  switch(
    measure,
    mse = apply((y_hat - y)^2, 3, mean),
    mae = apply(abs(y_hat - y), 3, mean)
  )
}

auc <- function(y, prob, weights = rep.int(1, nrow(y))) {
  if (is.matrix(y) || is.data.frame(y)) {
    ny <- nrow(y)
    auc(rep(c(0, 1), c(ny, ny)), c(prob, prob), as.vector(weights * y))
  } else {
    if (is.null(weights)) {
      rprob <- rank(prob)
      n1 <- sum(y)
      n0 <- length(y) - n1
      u <- sum(rprob[y == 1]) - n1 * (n1 + 1) / 2
      exp(log(u) - log(n1) - log(n0))
    } else {
      # randomize ties
      rprob <- stats::runif(length(prob))
      op <- order(prob, rprob)
      y <- y[op]
      weights <- weights[op]
      cw <- cumsum(weights)
      w1 <- weights[y == 1]
      cw1 <- cumsum(w1)
      wauc <- log(sum(w1 * (cw[y == 1] - cw1)))
      sumw1 <- cw1[length(cw1)]
      sumw2 <- cw[length(cw)] - sumw1
      exp(wauc - log(sumw1) - log(sumw2))
    }
  }
}
