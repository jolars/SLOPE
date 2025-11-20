#' Sorted L-One Penalized Estimation
#'
#' Fit a generalized linear model regularized with the
#' sorted L1 norm, which applies a
#' non-increasing regularization sequence to the
#' coefficient vector (\eqn{\beta}) after having sorted it
#' in decreasing order according  to its absolute values.
#'
#' `SLOPE()` solves the convex minimization problem
#' \deqn{
#'   f(\beta) + \alpha \sum_{i=j}^p \lambda_j |\beta|_{(j)},
#' }{
#'   f(\beta) + \alpha \sum \lambda_j |\beta|_(j),
#' }
#' where \eqn{f(\beta)} is a smooth and convex function and
#' the second part is the sorted L1-norm.
#' In ordinary least-squares regression,
#' \eqn{f(\beta)} is simply the squared norm of the least-squares residuals.
#' See section **Families** for specifics regarding the various types of
#' \eqn{f(\beta)} (model families) that are allowed in `SLOPE()`.
#'
#' By default, `SLOPE()` fits a path of models, each corresponding to
#' a separate regularization sequence, starting from
#' the null (intercept-only) model to an almost completely unregularized
#' model. These regularization sequences are parameterized using
#' \eqn{\lambda} and \eqn{\alpha}, with only \eqn{\alpha} varying along the
#' path. The length of the path can be manually, but will terminate
#' prematurely depending on
#' arguments `tol_dev_change`, `tol_dev_ratio`, and `max_variables`.
#' This means that unless these arguments are modified, the path is not
#' guaranteed to be of length `path_length`.
#'
#' @section Families:
#'
#' **Gaussian**
#'
#' The Gaussian model (Ordinary Least Squares) minimizes the following
#' objective:
#' \deqn{
#'   \frac{1}{2} \Vert y - X\beta\Vert_2^2
#' }{
#'   (1/(2))||y - X \beta||_2^2
#' }
#'
#' **Binomial**
#'
#' The binomial model (logistic regression) has the following objective:
#' \deqn{
#'   \sum_{i=1}^n \log\left(1+ \exp\left(
#'     - y_i \left(x_i^T\beta + \beta_0 \right) \right) \right)
#' }{
#'   \sum log(1+ exp(- y_i x_i^T \beta))
#' }
#' with \eqn{y \in \{-1, 1\}}{y in {-1, 1}}.
#'
#' **Poisson**
#'
#' In poisson regression, we use the following objective:
#'
#' \deqn{
#'   -\sum_{i=1}^n \left(y_i\left(
#'     x_i^T\beta + \beta_0\right) - \exp\left(x_i^T\beta + \beta_0
#'   \right)\right)
#' }{
#'   -\sum (y_i(x_i^T\beta + \beta_0) - exp(x_i^T\beta + \beta_0))
#' }
#'
#' **Multinomial**
#'
#' In multinomial regression, we minimize the full-rank objective
#' \deqn{
#'   -\sum_{i=1}^n\left(
#'     \sum_{k=1}^{m-1} y_{ik}(x_i^T\beta_k + \beta_{0,k})
#'     - \log\sum_{k=1}^{m-1} \exp\big(x_i^T\beta_k + \beta_{0,k}\big)
#'   \right)
#' }{
#'   -\sum(y_ik(x_i^T\beta_k + \beta_{0,k})
#'   - log(\sum exp(x_i^T\beta_k + \alpha_{0,k})))
#' }
#' with \eqn{y_{ik}} being the element in a \eqn{n} by \eqn{(m-1)} matrix, where
#' \eqn{m} is the number of classes in the response.
#'
#' @section Regularization Sequences:
#' There are multiple ways of specifying the `lambda` sequence
#' in `SLOPE()`. It is, first of all, possible to select the sequence manually
#' by
#' using a non-increasing
#' numeric vector, possibly of length one, as argument instead of a character.
#' The greater the differences are between
#' consecutive values along the sequence, the more clustering behavior
#' will the model exhibit. Note, also, that the scale of the \eqn{\lambda}
#' vector makes no difference if `alpha = NULL`, since `alpha` will be
#' selected automatically to ensure that the model is completely sparse at the
#' beginning and almost unregularized at the end. If, however, both
#' `alpha` and `lambda` are manually specified, then the scales of both do
#' matter, so make sure to choose them wisely.
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
#' for \eqn{i=1,\dots,p}, where \eqn{\Phi^{-1}}{\Phi^-1} is the quantile
#' function for the standard normal distribution and \eqn{q} is a parameter
#' that can be set by the user in the call to `SLOPE()`.
#'
#' **Gaussian**
#'
#' This penalty sequence is related to BH, such that
#' \deqn{
#'   \lambda_i = \lambda^{(\mathrm{BH})}_i
#'   \sqrt{1 + w(i-1)\cdot \mathrm{cumsum}(\lambda^2)_i},
#' }{
#'   \lambda_i = \lambda^(BH)_i \sqrt{1 + w(i-1) * cumsum(\lambda^2)_i},
#' }
#' for \eqn{i=1,\dots,p}, where \eqn{w(k) = 1/(n-k-1)}. We let
#' \eqn{\lambda_1 = \lambda^{(\mathrm{BH})}_1}{\lambda_1 = \lambda^(BH)_1} and
#' adjust the sequence to make sure that it's non-increasing.
#' Note that if \eqn{p} is large relative
#' to \eqn{n}, this option will result in a constant sequence, which is
#' usually not what you would want.
#'
#' **OSCAR**
#'
#' This sequence comes from Bondell and Reich and is a linear non-increasing
#' sequence, such that
#' \deqn{
#'   \lambda_i = \theta_1 + (p - i)\theta_2.
#' }
#' for \eqn{i = 1,\dots,p}. We use the parametrization from Zhong and Kwok
#' (2021) but use \eqn{\theta_1} and \eqn{\theta_2} instead of \eqn{\lambda_1}
#' and \eqn{\lambda_2} to avoid confusion and abuse of notation.
#'
#' **lasso**
#'
#' SLOPE is exactly equivalent to the
#' lasso when the sequence of regularization weights is constant, i.e.
#' \deqn{
#'   \lambda_i = 1
#' }
#' for \eqn{i = 1,\dots,p}. Here, again, we stress that the fact that
#' all \eqn{\lambda} are equal to one does not matter as long as
#' `alpha == NULL` since we scale the vector automatically.
#' Note that this option is only here for academic interest and
#' to highlight the fact that SLOPE is
#' a generalization of the lasso. There are more efficient packages, such as
#' **glmnet** and **biglasso**, for fitting the lasso.
#'
#' @section Solvers:
#'
#' There are currently three solvers available for SLOPE: Hybrid (Beck and
#' Teboulle 2009), proximal gradient descent (PGD), and FISTA (Beck and
#' Teboulle, 2009). The hybrid method is the preferred and generally
#' fastest method and is therefore the default for the Gaussian and
#' binomial families, but not currently available for multinomial and
#' disabled for Poisson due to convergence issues.
#'
#' @param x the design matrix, which can be either a dense
#'   matrix of the standard *matrix* class, or a sparse matrix
#'   inheriting from [Matrix::sparseMatrix]. Data frames will
#'   be converted to matrices internally.
#' @param y the response, which for `family = "gaussian"` must be numeric; for
#'   `family = "binomial"` or `family = "multinomial"`, it can be a factor.
#' @param family model family (objective); see **Families** for details.
#' @param intercept whether to fit an intercept
#' @param center whether to center predictors or not by their mean. Defaults
#'   to `TRUE` if `x` is dense and `FALSE` otherwise.
#' @param scale type of scaling to apply to predictors.
#'   - `"l1"` scales predictors to have L1 norms of one.
#'   - `"l2"` scales predictors to have L2 norms of one.#'
#'   - `"sd"` scales predictors to have a population standard deviation one.
#'   - `"none"` applies no scaling.
#' @param alpha scale for regularization path: either a decreasing numeric
#'   vector (possibly of length 1) or a character vector; in the latter case,
#'   the choices are:
#'   - `"path"`, which computes a regularization sequence
#'     where the first value corresponds to the intercept-only (null) model and
#'     the last to the almost-saturated model, and
#'   - `"estimate"`, which estimates a *single* `alpha`
#'     using Algorithm 5 in Bogdan et al. (2015).
#'
#'   When a value is manually entered for `alpha`, it will be scaled based
#'   on the type of standardization that is applied to `x`. For `scale = "l2"`,
#'   `alpha` will be scaled by \eqn{\sqrt n}. For `scale = "sd"` or `"none"`,
#'   alpha will be scaled by \eqn{n}, and for `scale = "l1"` no scaling is
#'   applied. Note, however, that the `alpha` that is returned in the
#'   resulting value is the **unstandardized** alpha.
#' @param path_length length of regularization path; note that the path
#'   returned may still be shorter due to the early termination criteria
#'   given by `tol_dev_change`, `tol_dev_ratio`, and `max_variables`.
#' @param lambda either a character vector indicating the method used
#'   to construct the lambda path or a numeric non-decreasing
#'   vector with length equal to the number
#'   of coefficients in the model; see section **Regularization sequences**
#'   for details.
#' @param alpha_min_ratio smallest value for `lambda` as a fraction of
#'   `lambda_max`; used in the selection of `alpha` when `alpha = "path"`.
#' @param q parameter controlling the shape of the lambda sequence, with
#'   usage varying depending on the type of path used and has no effect
#'   is a custom `lambda` sequence is used. Must be greater than `1e-6` and
#'   smaller than 1.
#' @param theta1 parameter controlling the shape of the lambda sequence
#'   when `lambda == "OSCAR"`. This parameter basically sets the intercept
#'   for the lambda sequence and is equivalent to \eqn{\lambda_1} in the
#'   original OSCAR formulation.
#' @param theta2 parameter controlling the shape of the lambda sequence
#'   when `lambda == "OSCAR"`. This parameter basically sets the slope
#'   for the lambda sequence and is equivalent to \eqn{\lambda_2} in the
#'   original OSCAR formulation.
#' @param max_passes maximum number of passes (outer iterations) for solver
#' @param diagnostics whether to save diagnostics from the solver
#'   (timings and other values depending on type of solver)
#' @param patterns whether to return the SLOPE pattern
#'   (cluster, ordering, and sign information) as a list of sparse
#'   matrices, one for each step on the path.
#' @param tol_dev_change the regularization path is stopped if the
#'   fractional change in deviance falls below this value; note that this is
#'   automatically set to 0 if a alpha is manually entered
#' @param tol_dev_ratio the regularization path is stopped if the
#'   deviance ratio \eqn{1 - \mathrm{deviance}/\mathrm{(null-deviance)}
#'   }{1 - deviance/(null deviance)} is above this threshold
#' @param max_variables criterion for stopping the path in terms of the
#'   maximum number of unique, nonzero coefficients in absolute value in model.
#'   For the multinomial family, this value will be multiplied internally with
#'   the number of levels of the response minus one.
#' @param tol stopping criterion for the solvers in terms of the relative
#'   duality gap
#' @param threads number of threads to use in the solver; if `NULL`, half
#'   of the available (logical) threads will be used
#' @param cd_type Type of coordinate descent to use, either `"cyclical"` or
#'  `"permuted"`. The former means that the cluster are cycled through
#'  in descending order of their coefficients' magnitudes, while the
#'  latter means that the clusters are permuted randomly for each pass.
#' @param verbosity DEPRECATED
#' @param prox_method DEPRECATED
#' @param screen DEPRECATED
#' @param screen_alg DEPRECATED
#' @param gamma relaxation mixing parameter, between 0 and 1. Has no effect if
#'   set to 0. If larger than 0, the solver will mix the coefficients
#'   from the ordinary SLOPE solutions with the coefficients from the
#'   relaxed solutions (fitting OLS on the SLOPE pattern).
#' @param tol_rel_gap DEPRECATED
#' @param tol_infeas DEPRECATED
#' @param tol_abs DEPRECATED
#' @param tol_rel relative DEPRECATED
#' @param tol_rel_coef_change DEPRECATED
#' @param solver type of solver use, either `"auto"`, `"hybrid"`, `"pgd"`, or
#'   `"fista"`; `"auto"` means that the solver is automatically selected,
#'   which currently means that `"hybrid"` is used for all objectives
#'   except multinomial ones, in which case FISTA (`"fista"`) is used.
#'
#' @return An object of class `"SLOPE"` with the following slots:
#' \item{coefficients}{
#'   a list of the coefficients from the
#'   model fit, not including the intercept.
#'   The coefficients are stored as sparse matrices.
#' }
#' \item{nonzeros}{
#'   a three-dimensional logical array indicating whether a
#'   coefficient was zero or not
#' }
#' \item{lambda}{
#'   the lambda vector that when multiplied by a value in `alpha`
#'   gives the penalty vector at that point along the regularization
#'   path
#' }
#' \item{alpha}{
#'   vector giving the (unstandardized) scaling of the lambda sequence
#' }
#' \item{class_names}{
#'   a character vector giving the names of the classes for binomial and
#'   multinomial families
#' }
#' \item{passes}{the number of passes the solver took at each step on the path}
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
#'   well as a measure of the infeasibility, time, and iteration; only
#'   available if `diagnostics = TRUE` in the call to [SLOPE()].
#' }
#' \item{n_observations}{
#'   the number of observations in the training data
#' }
#' \item{n_predictors}{
#'   the number of predictors in the training data
#' }
#' \item{call}{the call used for fitting the model}
#' @export
#'
#' @seealso [plot.SLOPE()], [plotDiagnostics()], [score()], [predict.SLOPE()],
#'   [trainSLOPE()], [coef.SLOPE()], [print.SLOPE()], [print.SLOPE()],
#'   [deviance.SLOPE()], [sortedL1Prox()]
#'
#' @references
#' Bogdan, M., van den Berg, E., Sabatti, C., Su, W., & Candès, E. J. (2015).
#' SLOPE -- adaptive variable selection via convex optimization. The Annals of
#' Applied Statistics, 9(3), 1103–1140.
#'
#' Larsson, J., Klopfenstein, Q., Massias, M., & Wallin, J. (2023). Coordinate
#' descent for SLOPE. In F. Ruiz, J. Dy, & J.-W. van de Meent (Eds.),
#' Proceedings of the 26th international conference on artificial
#' intelligence and statistics (Vol. 206, pp. 4802–4821). PMLR.
#' https://proceedings.mlr.press/v206/larsson23a.html
#'
#' Bondell, H. D., & Reich, B. J. (2008). Simultaneous Regression Shrinkage,
#' Variable Selection, and Supervised Clustering of Predictors with OSCAR.
#' Biometrics, 64(1), 115–123. JSTOR.
#'
#' Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2010).
#' Distributed Optimization and Statistical Learning via the Alternating
#' Direction Method of Multipliers. Foundations and Trends® in Machine Learning,
#' 3(1), 1–122.
#'
#' Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences,
#' 2(1), 183–202.
#'
#' @examples
#'
#' # Gaussian response, default lambda sequence
#' fit <- SLOPE(bodyfat$x, bodyfat$y)
#'
#' # Multinomial response, custom alpha and lambda
#' m <- length(unique(wine$y)) - 1
#' p <- ncol(wine$x)
#'
#' alpha <- 0.005
#' lambda <- exp(seq(log(2), log(1.8), length.out = p * m))
#'
#' fit <- SLOPE(
#'   wine$x,
#'   wine$y,
#'   family = "multinomial",
#'   lambda = lambda,
#'   alpha = alpha
#' )
SLOPE <- function(
  x,
  y,
  family = c("gaussian", "binomial", "multinomial", "poisson"),
  intercept = TRUE,
  center = c("mean", "min", "none"),
  scale = c("sd", "l1", "l2", "max_abs", "none"),
  alpha = c("path", "estimate"),
  lambda = c("bh", "gaussian", "oscar", "lasso"),
  alpha_min_ratio = if (NROW(x) < NCOL(x)) 1e-2 else 1e-4,
  path_length = 100,
  q = 0.1,
  theta1 = 1,
  theta2 = 0.5,
  tol_dev_change = 1e-5,
  tol_dev_ratio = 0.999,
  max_variables = NROW(x) + 1,
  solver = c("auto", "hybrid", "pgd", "fista", "admm"),
  max_passes = 1e6,
  tol = 1e-4,
  threads = 1,
  diagnostics = FALSE,
  patterns = FALSE,
  gamma = 1,
  cd_type = c("permuted", "cyclical"),
  tol_abs,
  tol_rel,
  tol_rel_gap,
  tol_infeas,
  tol_rel_coef_change,
  prox_method,
  screen,
  verbosity,
  screen_alg
) {
  if (!missing(tol_abs)) {
    warning("`tol_abs` argument is deprecated and has no effect")
  }

  if (!missing(tol_rel)) {
    warning("`tol_rel` argument is deprecated and has no effect")
  }

  if (!missing(tol_rel_gap)) {
    warning("`tol_rel_gap` argument is deprecated and has no effect")
  }

  if (!missing(tol_infeas)) {
    warning("`tol_infeas` argument is deprecated and has no effect")
  }

  if (!missing(tol_rel_coef_change)) {
    warning("`tol_rel_coef_change` argument is deprecated and has no effect")
  }

  if (!missing(screen)) {
    warning("`screen` argument is deprecated and has no effect")
  }

  if (!missing(screen_alg)) {
    warning("`screen_alg` argument is deprecated and has no effect")
  }

  if (!missing(verbosity)) {
    warning("`verbosity` argument is deprecated and has no effect")
  }

  if (!missing(prox_method)) {
    warning(
      "The 'method' argument is deprecated and ",
      "has no effect. It will be removed in a future version."
    )
  }

  ocall <- match.call()

  control <- processSlopeArgs(
    x,
    y,
    family,
    intercept,
    center,
    scale,
    alpha,
    lambda,
    alpha_min_ratio,
    path_length,
    q,
    theta1,
    theta2,
    tol_dev_change,
    tol_dev_ratio,
    max_variables,
    solver,
    max_passes,
    tol,
    diagnostics,
    patterns,
    threads,
    gamma,
    cd_type
  )

  x <- control$x
  y <- control$y
  alpha <- control$alpha

  is_sparse <- inherits(x, "sparseMatrix")
  is_big_matrix <- inherits(x, "big.matrix")

  fitSLOPE <- if (is_sparse) {
    sparseSLOPE
  } else if (is_big_matrix) {
    bigSLOPE
  } else {
    denseSLOPE
  }

  if (is_big_matrix) {
    x <- x@address
  }

  fit <- fitSLOPE(x, y, control)

  lambda <- fit$lambda
  alpha <- fit$alpha
  path_length <- length(alpha)
  intercepts <- fit$intercepts
  intercepts_scaled <- fit$intercepts_scaled
  beta <- fit$betas
  beta_scaled <- fit$betas_scaled
  # TODO: Do not return nonzeros; it's just wasting space.
  nonzeros <- lapply(beta, function(b) abs(b) > 0)
  coefficients <- beta

  names(coefficients) <- paste0("p", seq_along(beta))

  # check if maximum number of passes where reached anywhere
  passes <- fit$passes
  reached_max_passes <- passes >= max_passes

  if (any(reached_max_passes)) {
    reached_max_passes_where <- which(reached_max_passes)
    warning(
      "maximum number of passes reached at steps ",
      paste(reached_max_passes_where, collapse = ", "),
      "!"
    )
  }

  diagnostics <- if (diagnostics) setup_diagnostics(fit) else NULL

  patterns <- if (control$patterns) {
    fit$patterns
  } else {
    NULL
  }

  slope_class <- switch(
    control$family,
    gaussian = "GaussianSLOPE",
    binomial = "BinomialSLOPE",
    poisson = "PoissonSLOPE",
    multinomial = "MultinomialSLOPE"
  )

  structure(
    list(
      intercepts = intercepts,
      coefficients = coefficients,
      intercepts_scaled = intercepts_scaled,
      coefficients_scaled = beta_scaled,
      nonzeros = nonzeros,
      lambda = lambda,
      alpha = alpha[seq_along(beta)],
      variable_names = control$variable_names,
      class_names = control$class_names,
      passes = passes,
      deviance_ratio = drop(fit$deviance_ratio),
      null_deviance = fit$null_deviance,
      family = control$family,
      diagnostics = diagnostics,
      patterns = patterns,
      has_intercept = control$fit_intercept,
      n_observations = NROW(x),
      n_predictors = NCOL(x),
      call = ocall
    ),
    class = c(slope_class, "SLOPE")
  )
}

processSlopeArgs <- function(
  x,
  y,
  family = c("gaussian", "binomial", "multinomial", "poisson"),
  intercept = TRUE,
  center = c("mean", "min", "none"),
  scale = c("sd", "l1", "l2", "max_abs", "none"),
  alpha = c("path", "estimate"),
  lambda = c("bh", "gaussian", "oscar", "lasso"),
  alpha_min_ratio = if (NROW(x) < NCOL(x)) 1e-2 else 1e-4,
  path_length = 100,
  q = 0.1,
  theta1 = 1,
  theta2 = 0.5,
  tol_dev_change = 1e-5,
  tol_dev_ratio = 0.999,
  max_variables = NROW(x),
  solver = c("auto", "hybrid", "pgd", "fista", "admm"),
  max_passes = 1e6,
  tol = 1e-4,
  diagnostics = FALSE,
  patterns = TRUE,
  threads = NULL,
  gamma = 1,
  cd_type = c("permuted", "cyclical")
) {
  family <- match.arg(family)
  solver <- match.arg(solver)
  cd_type <- match.arg(cd_type)

  n <- NROW(x)
  p <- NCOL(x)

  stopifnot(
    is.null(alpha_min_ratio) ||
      (alpha_min_ratio > 0 && alpha_min_ratio < 1),
    max_passes > 0,
    q > 1e-6,
    q < 1,
    length(path_length) == 1,
    path_length >= 1,
    is.null(lambda) || is.character(lambda) || is.numeric(lambda),
    is.finite(max_passes),
    is.logical(diagnostics),
    is.logical(intercept),
    tol >= 0,
    is.finite(tol),
    theta1 >= 0,
    theta2 >= 0,
    is.finite(theta1),
    is.finite(theta2),
    gamma >= 0,
    gamma <= 1
  )

  if (solver == "admm") {
    warning("The ADMM solver has been removed; setting `solver` to `'auto'`")
    solver <- "auto"
  }

  centers <- double(0)
  scales <- double(0)

  if (is.character(center)) {
    center <- match.arg(center)
  } else if (is.logical(center) && length(center) == 1L) {
    center <- if (center) "mean" else "none"
  } else if (is.numeric(center)) {
    centers <- center
    center <- "manual"
  } else {
    stop(
      "`center` must be logical, a character, or a numeric of length ncol(x)"
    )
  }

  if (is.character(scale)) {
    scale <- match.arg(scale)
  } else if (is.logical(scale) && length(scale) == 1L) {
    scale <- if (scale) "sd" else "none"
  } else if (is.numeric(scale)) {
    scales <- scale
    scale <- "manual"
  } else {
    stop("`scale` must be logical, a character, or a numeric of length ncol(x)")
  }

  fit_intercept <- intercept

  # convert sparse x to dgCMatrix class from package Matrix.
  is_sparse <- inherits(x, "sparseMatrix")
  is_big_matrix <- inherits(x, "big.matrix")

  if (NROW(y) == 0) {
    stop("`y` is empty")
  }

  if (NROW(x) == 0) {
    stop("`x` is empty")
  }

  if (is_sparse) {
    x <- as_dgCMatrix(x)
  } else if (!is_big_matrix) {
    x <- as.matrix(x)
  }

  if (is.null(threads)) {
    threads <- -1
  } else {
    if (!is.numeric(threads) || length(threads) != 1L) {
      stop("`threads` must be a single numeric value")
    }
    threads <- as.integer(threads)
  }

  res <- preprocessResponse(family, y, fit_intercept)
  y <- as.matrix(res$y)
  class_names <- res$class_names
  m <- res$n_targets
  response_names <- res$response_names
  variable_names <- colnames(x)

  if (is.null(variable_names)) {
    variable_names <- paste0("V", seq_len(p))
  }
  if (is.null(response_names)) {
    response_names <- paste0("y", seq_len(m))
  }

  if (is.character(alpha)) {
    alpha <- match.arg(alpha)

    if (alpha == "path") {
      alpha_type <- "path"
      alpha <- double(0)
    } else if (alpha == "estimate") {
      if (family != "gaussian") {
        stop("`alpha = 'estimate'` can only be used if `family = 'gaussian'`")
      }

      alpha_type <- "estimate"
      alpha <- 0
      path_length <- 1
    }
  } else {
    alpha <- as.double(alpha)
    alpha_type <- "path"

    alpha <- as.double(alpha)
    path_length <- length(alpha)
  }

  if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)

    if (lambda_type == "bhq") {
      warning(
        "'bhq' option to argument lambda has been depracted and will",
        "will be defunct in the next release; please use 'bh' instead"
      )
    }

    lambda <- double(0)
  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)
  }

  control <- list(
    x = x,
    y = as.matrix(y),
    alpha = alpha,
    alpha_min_ratio = alpha_min_ratio,
    alpha_type = alpha_type,
    centering_type = center,
    centers = centers,
    class_names = class_names,
    diagnostics = diagnostics,
    patterns = patterns,
    family = family,
    fit_intercept = fit_intercept,
    gamma = gamma,
    lambda = lambda,
    lambda_type = lambda_type,
    max_passes = max_passes,
    max_variables = max_variables,
    path_length = path_length,
    q = q,
    response_names = response_names,
    scales = scales,
    scaling_type = scale,
    solver = solver,
    theta1 = theta1,
    theta2 = theta2,
    tol = tol,
    tol_dev_change = tol_dev_change,
    tol_dev_ratio = tol_dev_ratio,
    variable_names = variable_names,
    threads = threads,
    cd_type = cd_type
  )

  control
}
