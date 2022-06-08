#' Model objects for model tuning with caret (deprecated)
#'
#' This function can be used in a call to [caret::train()] to enable
#' model tuning using caret. Note that this function does not properly work
#' with sparse feature matrices and standardization due to the way
#' resampling is implemented in caret. So for these cases, please
#' check out [trainSLOPE()] instead.
#'
#' @return A model description list to be used in the `method` argument
#' in [caret::train()].
#'
#' @seealso [caret::train()], [trainSLOPE()], [SLOPE()]
#' @family model-tuning
#'
#' @export
caretSLOPE <- function() {
  .Deprecated()

  list(
    label = "SLOPE",
    library = c("SLOPE", "Matrix"),
    type = c("Regression", "Classification"),
    parameters = data.frame(
      parameter = c("alpha", "q"),
      class = rep("numeric", 2),
      label = c("alpha", "False discovery rate")
    ),
    grid = function(x, y, len = NULL, search = "grid") {
      s <- 0.2 * (1 - 1 / len)
      q <- seq(-s, s, length.out = len) + 0.2

      numLev <- if (is.character(y) | is.factor(y)) length(levels(y)) else NA

      if (!is.na(numLev)) {
        fam <- ifelse(numLev > 2, "multinomial", "binomial")
      } else {
        fam <- "gaussian"
      }

      fit <- SLOPE::SLOPE(x,
        y,
        family = fam,
        path_length = len + 2,
        scale = FALSE,
        center = FALSE
      )

      alpha <- fit$alpha
      alpha <- alpha[-c(1, length(alpha))]
      alpha <- alpha[1:min(length(alpha), len)]

      if (search == "grid") {
        out <- expand.grid(
          alpha = alpha,
          q = q
        )
      } else {
        q <- stats::runif(0.01, 0.4, len)

        alpha <- exp(stats::runif(
          log(min(alpha)),
          log(max(alpha)),
          len
        ))

        out <- data.frame(
          alpha = alpha,
          q = q
        )
      }

      out
    },
    loop = function(grid) {
      q <- unique(grid$q)
      loop <- data.frame(q = q)
      loop$alpha <- NA

      submodels <- vector("list", length = length(q))

      for (i in seq_along(q)) {
        np <- grid[grid$q == q[i], "alpha"]
        loop$alpha[loop$q == q[i]] <- np[which.max(np)]
        submodels[[i]] <- data.frame(alpha = np[-which.max(np)])
      }

      list(loop = loop, submodels = submodels)
    },
    fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
      dots <- list(...)

      numLev <- if (is.character(y) | is.factor(y)) {
        length(levels(y))
      } else {
        NA
      }

      if (all(names(dots) != "family")) {
        if (!is.na(numLev)) {
          fam <- ifelse(numLev > 2, "multinomial", "binomial")
        } else {
          fam <- "gaussian"
        }

        dots$family <- fam
      }

      dots$x <- x
      dots$y <- y
      dots$q <- param$q
      dots$scale <- FALSE
      dots$center <- FALSE
      dots$tol_dev_change <- 0
      dots$tol_dev_ratio <- 1

      out <- do.call(SLOPE::SLOPE, dots)

      if (!is.na(param$alpha[1])) {
        out$alphaOpt <- param$alpha[1]
      }

      out
    },
    predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
      library(SLOPE)

      if (!is.matrix(newdata)) {
        newdata <- Matrix::as.matrix(newdata)
      }

      if (length(modelFit$obsLevels) < 2) {
        out <- stats::predict(modelFit,
          newdata,
          alpha = modelFit$alphaOpt
        )
      } else {
        out <- stats::predict(modelFit,
          newdata,
          alpha = modelFit$alphaOpt,
          type = "class"
        )
      }

      if (is.matrix(out)) {
        out <- out[, 1]
      }

      if (!is.null(submodels)) {
        if (length(modelFit$obsLevels) < 2) {
          tmp <- stats::predict(modelFit,
            newdata,
            alpha = submodels$alpha
          )
          tmp <- as.list(as.data.frame(tmp))
        } else {
          tmp <- stats::predict(modelFit,
            newdata,
            alpha = submodels$alpha,
            type = "class"
          )
          tmp <- if (is.matrix(tmp)) {
            as.data.frame(tmp, stringsAsFactors = FALSE)
          } else {
            as.character(tmp)
          }
          tmp <- as.list(tmp)
        }
        out <- c(list(out), tmp)
      }
      out
    },
    prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
      obsLevels <- if ("class_names" %in% names(modelFit)) {
        modelFit$class_names
      } else {
        NULL
      }

      probs <- stats::predict(modelFit,
        Matrix::as.matrix(newdata),
        alpha = modelFit$alpha,
        type = "response"
      )

      if (length(obsLevels) == 2) {
        probs <- as.vector(probs)
        probs <- as.data.frame(cbind(1 - probs, probs))
        colnames(probs) <- modelFit$obsLevels
      } else {
        probs <- as.data.frame(probs[, , 1, drop = FALSE])
        names(probs) <- modelFit$obsLevels
      }

      if (!is.null(submodels)) {
        tmp <- stats::predict(modelFit,
          Matrix::as.matrix(newdata),
          alpha = submodels$alpha,
          type = "response"
        )

        if (length(obsLevels) == 2) {
          tmp <- as.list(as.data.frame(tmp))
          tmp <- lapply(tmp,
            function(x, lev) {
              x <- as.vector(x)
              tmp <- data.frame(1 - x, x)
              names(tmp) <- lev
              tmp
            },
            lev = modelFit$obsLevels
          )
        } else {
          tmp <- apply(tmp, 3, function(x) data.frame(x))
        }
        probs <- if (is.list(tmp)) c(list(probs), tmp) else list(probs, tmp)
      }
      probs
    },
    predictors = function(x, alpha = NULL, ...) {
      if (is.null(alpha)) {
        if (length(alpha) > 1) {
          stop("Only one value of alpha is allowed right now")
        }

        if (!is.null(x$alphaOpt)) {
          alpha <- x$alphaOpt
        } else {
          stop("must supply a value of alpha")
        }
      }

      allVar <- rownames(stats::coef(x, simplify = FALSE))

      out <- apply(abs(allVar) > 0, 3, sum)

      out <- unique(out)
      if (length(out) > 0) {
        out <- out[!is.na(out)]
        out <- allVar[out]
      }

      out
    },
    varImp = function(object, alpha = NULL, ...) {
      if (is.null(alpha)) {
        if (length(alpha) > 1) {
          stop("Only one value of alpha is allowed right now")
        }

        if (!is.null(object$alphaOpt)) {
          alpha <- object$alphaOpt
        } else {
          stop("must supply a value of alpha")
        }
      }

      beta <- stats::coef(object, alpha = alpha, simplify = TRUE)
      beta <- as.data.frame(beta)
      out <- as.data.frame(Overall = beta[, 1])
      out <- abs(out[rownames(out) != "(Intercept)", , drop = FALSE])
      out
    },
    levels = function(x) {
      if (any(names(x) == "obsLevels")) {
        x$obsLevels
      } else {
        NULL
      }
    },
    sort = function(x) {
      x[order(-x$alpha, x$q), ]
    },
    trim = function(x) {
      x$call <- NULL
      x
    },
    tags = c(
      "Generalized Linear Model",
      "Implicit Feature Selection",
      "L1 Regularization",
      "Linear Classifier",
      "Linear Regression"
    )
  )
}
