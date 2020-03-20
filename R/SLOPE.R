#' Generalized linear models regularized with the Sorted L1 Norm
#'
#' Fit a generalized linear model regularized with the
#' SLOPE (Sorted L-One Penalized Estimation) norm, which applies a
#' decreasing \eqn{\lambda} (penalty sequence) to the
#' coefficient vector (\eqn{\beta}) after having sorted it
#' in decreasing order according  to its absolute values.
#'
#' `SLOPE()` tries to minimize the following composite objective, given
#' in Lagrangian form.
#' \deqn{
#'   f(\beta) + \sigma \sum_{i=j}^p \lambda_j |\beta|_{(j)},
#' }{
#'   f(\beta) + \sigma \sum \lambda |\beta|(j),
#' }
#' where \eqn{f(\beta)} is a smooth, convex function of \eqn{\beta}, whereas
#' the second part is the SLOPE norm, which is convex but non-smooth.
#' In ordinary least-squares regression, for instance,
#' \eqn{f(\beta)} is simply the squared norm of the least-squares residuals.
#' See section **Families** for specifics regarding the various types of
#' \eqn{f(\beta)} (model families) that are allowed in `SLOPE()`.
#'
#' By default, `SLOPE()` fits a path of `lambda` sequences, starting from
#' the null (intercept-only) model to an almost completely unregularized
#' model. The path will end prematurely, however, if the criteria
#' related to *any* of the
#' arguments `tol_dev_change`, `tol_dev_ratio`, or `max_variables`
#' are reached before the path is complete. This means that unless these
#' arguments are modified, the path is not guaranteed to be of
#' length `n_sigma`.
#'
#' @section Families:
#'
#' **Gaussian**
#'
#' The Gaussian model (Ordinary Least Squares) minimizes the following
#' objective.
#' \deqn{
#'   ||y - X\beta||_2^2
#' }
#'
#' **Binomial**
#'
#' The binomial model (logistic regression) has the following objective.
#' \deqn{
#'   \sum_{i=1}^n \log\left(1+ \exp\left(- y_i \left(x_i^T\beta + \alpha \right) \right) \right)
#' }{
#'   \sum log(1+ exp(- y_i x_i^T \beta))
#' }
#' with \eqn{y \in \{-1, 1\}}{y in {-1, 1}}.
#'
#' **Poisson**
#'
#' In poisson regression, we use the following objective.
#'
#' \deqn{
#'   -\sum_{i=1}^n \left(y_i\left(x_i^T\beta + \alpha\right) - \exp\left(x_i^T\beta + \alpha\right)\right)
#' }{
#'   -\sum (y_i(x_i^T\beta + \alpha) - exp(x_i^T\beta + \alpha))
#' }
#'
#' **Multinomial**
#'
#' In multinomial regression, we minimize the full-rank objective
#' \deqn{
#'   -\sum_{i=1}^n\left( \sum_{k=1}^{m-1} y_{ik}(x_i^T\beta_k + \alpha_k)
#'                      - \log\sum_{k=1}^{m-1} \exp(x_i^T\beta_k + \alpha_k) \right)
#' }{
#'   -\sum(y_ik(x_i^T\beta_k + \alpha_k) - log(\sum exp(x_i^T\beta_k + \alpha_k)))
#' }
#' with \eqn{y_{ik}} being the element in a \eqn{n} by \eqn{(m-1)} matrix, where
#' \eqn{m} is the number of classes in the response.
#'
#' @section Regularization sequences:
#' There are multiple ways of specifying the `lambda` sequence
#' in `SLOPE()`. It is, first of all, possible to select the sequence manually by
#' using a non-increasing
#' numeric vector as argument instead of a character.
#' If all `lambda` are the same value, this will
#' lead to the ordinary lasso penalty. The greater the differences are between
#' consecutive values along the sequence, the more clustering behavior
#' will the model exhibit. Note, also, that the scale of the \eqn{\lambda}
#' vector makes no difference if `sigma = NULL`, since `sigma` will be
#' selected automatically to ensure that the model is completely sparse at the
#' beginning and almost unregularized at the end. If, however, both
#' `sigma` and `lambda` are manually specified, both of the scales will
#' matter.
#'
#' Instead of choosing the sequence manually, one of the following
#' automatically generated sequences may be chosen.
#'
#' **BH (Benjamini--Hochberg)**
#'
#' If `lambda = "bh"`, the sequence used is that referred to
#' as \eqn{\lambda^{(\mathrm{BH})}}{\lambda^(BH)} by Bogdan et al, which sets
#' \eqn{\lambda} according to
#' \deqn{
#'   \lambda_i = \Phi^{-1}(1 - iq/(2p)),
#' }{
#'   \lambda_i = \Phi^-1(1 - iq/(2p)),
#' }
#' where \eqn{\Phi^{-1}}{\Phi^-1} is the quantile function for the standard
#' normal distribution and \eqn{q} is a parameter that can be
#' set by the user in the call to `SLOPE()`.
#'
#' **Gaussian**
#'
#' This penalty sequence is related to BH, such that
#' \deqn{
#'   \lambda_i = \lambda^{(\mathrm{BH})}_i \sqrt{1 + w(i-1)\cdot \mathrm{cumsum}(\lambda^2)_i},
#' }{
#'   \lambda_i = \lambda^(BH)_i \sqrt{1 + w(i-1) * cumsum(\lambda^2)_i},
#' }
#' where \eqn{w(k) = 1/(n-k-1)}. We let
#' \eqn{\lambda_1 = \lambda^{(\mathrm{BH})}_1}{\lambda_1 = \lambda^(BH)_1} and
#' adjust the sequence to make sure that it's non-increasing.
#' Note that if \eqn{p} is large relative
#' to \eqn{n}, this option will result in a constant sequence, which is
#' usually not what you would want.
#'
#' **OSCAR**
#'
#' This sequence comes from Bondell and Reich and is a linearly non-increasing
#' sequence such that
#' \deqn{
#'   \lambda_i = q(p - i) + 1.
#' }
#'
#' @param x the feature matrix, which can be either a dense
#'   matrix of the standard *matrix* class, or a sparse matrix
#'   inheriting from [Matrix::sparseMatrix] Data frames will
#'   be converted to matrices internally.
#' @param y the response. For Gaussian models this must be numeric; for
#'   binomial models, it can be a factor.
#' @param family response type. See **Families** for details.
#' @param intercept whether to fit an intercept
#' @param center whether to center predictors or not by their mean. Defaults
#'   to true if dense matrix, false otherwise.
#' @param scale type of scaling to apply to predictors, `"l1"` scales
#'   predictors to have L1-norm of one, `"l2"` scales predictors to have
#'   L2-norm one, `"sd"` scales predictors to have standard deviation one.
#' @param sigma scale of lambda sequence
#' @param n_sigma length of regularization path
#' @param lambda either a character vector indicating the method used
#'   to construct the lambda path or a vector with length equal to the number
#'   of coefficients in the model
#' @param lambda_min_ratio smallest value for `lambda` as a fraction of
#'   `lambda_max`
#' @param q shape of lambda sequence
#' @param max_passes maximum number of passes for optimizer
#' @param diagnostics should diagnostics be saved for the model fit (timings,
#'   primal and dual objectives, and infeasibility)
#' @param screening whether the strong rule for SLOPE be used to screen
#'   variables for inclusion
#' @param verbosity level of verbosity for displaying output from the
#'   program. Setting this to 1 displays basic information on the path level,
#'   2 a little bit more information on the path level, and 3 displays
#'   information from the solver.
#' @param tol_dev_change the regularization path is stopped if the
#'   fractional change in deviance falls below this value. Note that this is
#'   automatically set to 0 if a sigma is manually entered
#' @param tol_dev_ratio the regularization path is stopped if the
#'   deviance ratio
#'   \eqn{1 - \mathrm{deviance}/\mathrm{(null-deviance)}}{1 - deviance/(null deviance)}
#'   is above this threshold
#' @param max_variables criterion for stopping the path in terms of the
#'   maximum number of unique, nonzero
#'   coefficients in absolute value in model
#' @param tol_rel_gap stopping criterion for the duality gap
#' @param tol_infeas stopping criterion for the level of infeasibility
#' @param tol_abs absolute tolerance criterion for ADMM solver (used for
#'   Gaussian dense designs)
#' @param tol_rel relative tolerance criterion for ADMM solver (used for
#'   Gaussian dense designs)
#' @param X deprecated. please use `x` instead
#' @param fdr deprecated. please use `q` instead
#' @param normalize deprecated. please use `scale` and `center` instead
#' @param solver deprecated
#'
#' @return An object of class `"SLOPE"` with the following slots:
#' \item{coefficients}{
#'   a three-dimensional array of the coefficients from the
#'   model fit, including the intercept if it was fit.
#'   There is one row for each coefficient, one column
#'   for each target (dependent variable), and
#'   one slice for each penalty.
#' }
#' \item{nonzeros}{
#'   a three-dimensional boolean array indicating whether a
#'   coefficient was zero or not
#' }
#' \item{lambda}{
#'   the lambda vector that when multiplied by a value in `sigma`
#'   gives the penalty vector at that point along the regularization
#'   path
#' }
#' \item{sigma}{the vector of sigma, indicating the scale of the lambda vector}
#' \item{class_names}{
#'   a character vector giving the names of the classes for binomial and
#'   multinomial families
#' }
#' \item{passes}{the number of passes the solver took at each path}
#' \item{violations}{the number of violations of the screening rule}
#' \item{active_sets}{
#'   a list where each element indicates the indices of the
#'   coefficients that were active at that point in the regularization path
#' }
#' \item{unique}{
#'   the number of unique predictors (in absolute value)
#' }
#' \item{deviance_ratio}{
#'   the deviance ratio (as a fraction of 1)
#' }
#' \item{null_deviance}{
#'   the deviance of the null (intercept-only) model
#' }
#' \item{family}{
#'   the name of the family used in the model fit
#' }
#' \item{diagnostics}{
#'   a `data.frame` of objective values for the primal and dual problems, as
#'   well as a measure of the infeasibility, time, and iteration. Only
#'   available if `diagnostics = TRUE` in the call to [SLOPE()].
#' }
#' \item{call}{the call used for fitting the model}
#' @export
#'
#' @seealso [plot.SLOPE()], [plotDiagnostics()], [score()], [predict.SLOPE()],
#'   [trainSLOPE()], [coef.SLOPE()], [print.SLOPE()]
#'
#' @references
#' Bogdan, M., van den Berg, E., Sabatti, C., Su, W., & Candès, E. J. (2015).
#' SLOPE -- adaptive variable selection via convex optimization. The Annals of
#' Applied Statistics, 9(3), 1103–1140. <https://doi.org/10/gfgwzt>
#'
#' Bondell, H. D., & Reich, B. J. (2008). Simultaneous Regression Shrinkage,
#' Variable Selection, and Supervised Clustering of Predictors with OSCAR.
#' Biometrics, 64(1), 115–123. JSTOR.
#' <https://doi.org/10.1111/j.1541-0420.2007.00843.x>
#'
#'
#' @examples
#'
#' # Gaussian response, default lambda sequence
#'
#' fit <- SLOPE(bodyfat$x, bodyfat$y)
#'
#' # Binomial response, BH-type lambda sequence
#'
#' fit <- SLOPE(heart$x, heart$y, family = "binomial", lambda = "bh")
#'
#' # Poisson response, OSCAR-type lambda sequence
#'
#' fit <- SLOPE(abalone$x,
#'              abalone$y,
#'              family = "poisson",
#'              lambda = "oscar",
#'              q = 0.4)
#'
#' # Multinomial response, custom sigma and lambda
#'
#' m <- length(unique(wine$y)) - 1
#' p <- ncol(wine$x)
#'
#' sigma <- 0.005
#' lambda <- exp(seq(log(2), log(1.8), length.out = p*m))
#'
#' fit <- SLOPE(wine$x,
#'              wine$y,
#'              family = "multinomial",
#'              lambda = lambda,
#'              sigma = sigma)
#'
SLOPE <- function(x,
                  y,
                  family = c("gaussian", "binomial", "multinomial", "poisson"),
                  intercept = TRUE,
                  center = !inherits(x, "sparseMatrix"),
                  scale = c("l2", "l1", "sd", "none"),
                  sigma = c("path", "estimate"),
                  lambda = c("gaussian", "bh", "oscar"),
                  lambda_min_ratio = if (n < p) 1e-2 else 1e-4,
                  n_sigma = 100,
                  q = 0.1*min(1, n/p),
                  screening = TRUE,
                  tol_dev_change = 1e-5,
                  tol_dev_ratio = 0.995,
                  tol_abs = 1e-5,
                  tol_rel = 1e-4,
                  max_variables = n*m,
                  max_passes = 1e6,
                  tol_rel_gap = 1e-5,
                  tol_infeas = 1e-3,
                  diagnostics = FALSE,
                  verbosity = 0,
                  X,
                  fdr,
                  normalize,
                  solver
) {

  if (!missing(X)) {
    x <- X
    warning("'X' argument is deprecated; please use 'x' instead")
  }

  if (!missing(fdr)) {
    warning("'fdr' argument is deprecated; please use 'q' instead")
    q <- fdr
  }

  if (!missing(normalize)) {
    warning("'normalize' argument is deprecated; please use 'scale' and",
            "'center' instead. 'scale' has been set to 'l2', and",
            "'center' to TRUE")
    center <- TRUE
    scale <- "l2"
  }

  if (!missing(solver)) {
    warning("'solver' argument is deprecated")
  }

  ocall <- match.call()

  family <- match.arg(family)

  if (is.character(sigma)) {
    sigma <- match.arg(sigma)
    if (sigma == "estimate" && family != "gaussian")
      stop("sigma == 'estimate' can only be used if family == 'gaussian'")
  } else {
    sigma <- as.double(sigma)
  }

  if (is.character(scale)) {
    scale <- match.arg(scale)
  } else if (is.logical(scale) && length(scale) == 1L) {
    scale <- ifelse(scale, "l2", "none")
  } else {
    stop("'scale' must be logical or a character")
  }

  n <- NROW(x)
  p <- NCOL(x)

  stopifnot(
    is.null(lambda_min_ratio) ||
      (lambda_min_ratio > 0 && lambda_min_ratio < 1),
    max_passes > 0,
    q > 0,
    q < 1,
    length(n_sigma) == 1,
    n_sigma >= 1,
    is.null(lambda) || is.character(lambda) || is.numeric(lambda),
    is.finite(max_passes),
    is.logical(diagnostics),
    is.logical(intercept),
    tol_rel_gap >= 0,
    tol_infeas >= 0,
    tol_abs >= 0,
    tol_rel >= 0,
    is.logical(center)
  )

  fit_intercept <- intercept

  # convert sparse x to dgCMatrix class from package Matrix.
  is_sparse <- inherits(x, "sparseMatrix")

  if (NROW(y) != NROW(x))
    stop("the number of samples in 'x' and 'y' must match")

  if (NROW(y) == 0)
    stop("y is empty")

  if (NROW(x) == 0)
    stop("x is empty")

  if (anyNA(y) || anyNA(x))
    stop("missing values are not allowed")

  if (is_sparse) {
    x <- methods::as(x, "dgCMatrix")
  } else {
    x <- as.matrix(x)
  }

  if (is_sparse && center)
    stop("centering would destroy sparsity in x (predictors)")

  res <- preprocessResponse(family, y)
  y <- as.matrix(res$y)
  y_center <- res$y_center
  y_scale <- res$y_scale
  class_names <- res$class_names
  m <- n_targets <- res$n_targets
  response_names <- res$response_names
  variable_names <- colnames(x)

  if (is.null(variable_names))
    variable_names <- paste0("V", seq_len(p))
  if (is.null(response_names))
    response_names <- paste0("y", seq_len(m))

  if (sigma == "path") {
    sigma_type <- "auto"
    sigma <- double(n_sigma)
  } else if (sigma == "estimate") {
    sigma_type <- "estimate"
    sigma <- NULL
    if (n_sigma > 1)
      warning("'n_sigma' ignored since sigma == 'estimate'")
  } else {
    sigma_type <- "user"

    sigma <- as.double(sigma)
    n_sigma <- length(sigma)

    stopifnot(n_sigma > 0)

    # do not stop path early if user requests specific sigma
    tol_dev_change <- 0
    tol_dev_ratio <- 1
    max_variables <- (NCOL(x) + intercept)*m
  }

  n_lambda <- m*p

  if (is.null(lambda)) {
    lambda_type <- "bh"
    lambda <- double(n_lambda)
  } else if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)
    lambda <- double(n_lambda)
  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)

    if (length(lambda) != n_lambda)
      stop("lambda sequence must be as long as there are variables")

    if (is.unsorted(rev(lambda)))
      stop("lambda sequence must be non-increasing")

    if (any(lambda < 0))
      stop("lambda sequence cannot contain negative values")
  }

  control <- list(family = family,
                  fit_intercept = fit_intercept,
                  is_sparse = is_sparse,
                  scale = scale,
                  center = center,
                  n_sigma = n_sigma,
                  n_targets = n_targets,
                  screening = screening,
                  sigma = sigma,
                  sigma_type = sigma_type,
                  lambda = lambda,
                  lambda_type = lambda_type,
                  lambda_min_ratio = lambda_min_ratio,
                  q = q,
                  y_center = y_center,
                  y_scale = y_scale,
                  max_passes = max_passes,
                  diagnostics = diagnostics,
                  verbosity = verbosity,
                  max_variables = max_variables,
                  tol_dev_change = tol_dev_change,
                  tol_dev_ratio = tol_dev_ratio,
                  tol_rel_gap = tol_rel_gap,
                  tol_infeas = tol_infeas,
                  tol_abs = tol_abs,
                  tol_rel = tol_rel)

  fitSLOPE <- if (is_sparse) sparseSLOPE else denseSLOPE

  if (sigma_type %in% c("path", "user")) {
    if (intercept) {
      fit <- fitSLOPE(cbind(1, x), y, control)
    } else {
      fit <- fitSLOPE(x, y, control)
    }
  } else {
    # Estimate the noise level, if possible.
    if (is.null(sigma) && n >= p + 30)
      sigma <- estimateNoise(x, y)

    # Run the solver, iteratively if necessary.
    if (is.null(sigma)) {
      # Run Algorithm 5 of Section 3.2.3.
      selected <- integer(0)
      repeat {
        selected_prev <- selected
        sigma <- estimateNoise(x[, selected, drop = FALSE], y, intercept)
        control$sigma <- sigma

        if (intercept) {
          fit <- fitSLOPE(cbind(1, x), y, control)
        } else {
          fit <- fitSLOPE(x, y, control)
        }
        # result <- SLOPE_solver_call(solver, X, y, tail(sigma, 1) * lambda)
        selected <- which(abs(drop(fit$betas)) > 0)

        if (fit_intercept)
          selected <- selected[-1, , , drop = FALSE]

        if (identical(selected, selected_prev))
          break

        if (length(selected) + 1 >= n)
          stop("selected >= n-1 variables; cannot estimate variance")
      }
    } else {
      control$sigma <- sigma
      if (intercept) {
        fit <- fitSLOPE(cbind(1, x), y, control)
      } else {
        fit <- fitSLOPE(x, y, control)
      }
    }
  }

  lambda <- fit$lambda
  sigma <- fit$sigma
  n_sigma <- length(sigma)
  active_sets <- lapply(drop(fit$active_sets), function(x) drop(x) + 1)
  beta <- fit$betas
  nonzeros <- apply(beta, c(2, 3), function(x) abs(x) > 0)
  coefficients <- beta

  if (fit_intercept) {
    nonzeros <- nonzeros[-1, , , drop = FALSE]
    dimnames(coefficients) <- list(c("(Intercept)", variable_names),
                                   response_names[1:n_targets],
                                   paste0("p", seq_len(n_sigma)))
  } else {
    dimnames(coefficients) <- list(variable_names,
                                   response_names[1:n_targets],
                                   paste0("p", seq_len(n_sigma)))
  }

  diagnostics <- if (diagnostics) setupDiagnostics(fit) else NULL

  slope_class <- switch(family,
                        gaussian = "GaussianSLOPE",
                        binomial = "BinomialSLOPE",
                        poisson = "PoissonSLOPE",
                        multinomial = "MultinomialSLOPE")

  structure(list(coefficients = coefficients,
                 nonzeros = nonzeros,
                 lambda = lambda,
                 sigma = sigma,
                 class_names = class_names,
                 passes = fit$passes,
                 violations = fit$violations,
                 active_sets = active_sets,
                 unique = fit$n_unique,
                 deviance_ratio = as.vector(fit$deviance_ratio),
                 null_deviance = fit$null_deviance,
                 family = family,
                 diagnostics = diagnostics,
                 call = ocall),
            class = c(slope_class, "SLOPE"))
}
