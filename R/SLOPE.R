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
#' There are currently two solvers available for SLOPE: FISTA (Beck and
#' Teboulle 2009) and ADMM (Boyd et al. 2008). FISTA is available for
#' families but ADMM is currently only available for `family = "gaussian"`.
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
#' @param prox_method method for calculating the proximal operator for
#'   the Sorted L1 Norm (the SLOPE penalty). Please see [sortedL1Prox()] for
#'   more information.
#' @param max_passes maximum number of passes (outer iterations) for solver
#' @param diagnostics whether to save diagnostics from the solver
#'   (timings and other values depending on type of solver)
#' @param screen whether to use predictor screening rules (rules that allow
#'   some predictors to be discarded prior to fitting), which improve speed
#'   greatly when the number of predictors is larger than the number
#'   of observations.
#' @param screen_alg what type of screening algorithm to use.
#'   - `"strong"` uses the set from the strong screening rule and check
#'     against the full set
#'   - `"previous"` first fits with the previous active set, then checks
#'     against the strong set, and finally against the full set if there are
#'     no violations in the strong set
#' @param verbosity level of verbosity for displaying output from the
#'   program. Setting this to 1 displays basic information on the path level,
#'   2 a little bit more information on the path level, and 3 displays
#'   information from the solver.
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
#' @param tol_rel_gap stopping criterion for the duality gap; used only with
#'   FISTA solver.
#' @param tol_infeas stopping criterion for the level of infeasibility; used
#'   with FISTA solver and KKT checks in screening algorithm.
#' @param tol_abs absolute tolerance criterion for ADMM solver
#' @param tol_rel relative tolerance criterion for ADMM solver
#' @param tol_rel_coef_change relative tolerance criterion for change
#'   in coefficients between iterations, which is reached when
#'   the maximum absolute change in any coefficient divided by the maximum
#'   absolute coefficient size is less than this value.
#' @param sigma deprecated; please use `alpha` instead
#' @param n_sigma deprecated; please use `path_length` instead
#' @param lambda_min_ratio deprecated; please use `alpha_min_ratio` instead
#' @param solver type of solver use, either `"fista"` or `"admm"`;
#'   all families currently support FISTA but only `family = "gaussian"`
#'   supports ADMM.
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
#' \item{violations}{
#'   the number of violations of the screening rule at each step on the path;
#'   only available if `diagnostics = TRUE` in the call to [SLOPE()].
#' }
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
#'   well as a measure of the infeasibility, time, and iteration; only
#'   available if `diagnostics = TRUE` in the call to [SLOPE()].
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
#' Applied Statistics, 9(3), 1103–1140. \doi{10/gfgwzt}
#'
#' Bondell, H. D., & Reich, B. J. (2008). Simultaneous Regression Shrinkage,
#' Variable Selection, and Supervised Clustering of Predictors with OSCAR.
#' Biometrics, 64(1), 115–123. JSTOR.
#' \doi{10.1111/j.1541-0420.2007.00843.x}
#'
#' Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2010).
#' Distributed Optimization and Statistical Learning via the Alternating
#' Direction Method of Multipliers. Foundations and Trends® in Machine Learning,
#' 3(1), 1–122. \doi{10.1561/2200000016}
#'
#' Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding
#' Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences,
#' 2(1), 183–202. \doi{10.1137/080716542}
#'
#' @examples
#'
#' # Gaussian response, default lambda sequence
#' fit <- SLOPE(bodyfat$x, bodyfat$y)
#'
#' # Poisson response, OSCAR-type lambda sequence
#' fit <- SLOPE(
#'   abalone$x,
#'   abalone$y,
#'   family = "poisson",
#'   lambda = "oscar",
#'   theta1 = 1,
#'   theta2 = 0.9
#' )
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
SLOPE <- function(x,
                  y,
                  family = c("gaussian", "binomial", "multinomial", "poisson"),
                  intercept = TRUE,
                  center = !inherits(x, "sparseMatrix"),
                  scale = c("l2", "l1", "sd", "none"),
                  alpha = c("path", "estimate"),
                  lambda = c("bh", "gaussian", "oscar", "lasso"),
                  alpha_min_ratio = if (NROW(x) < NCOL(x)) 1e-2 else 1e-4,
                  path_length = if (alpha[1] == "estimate") 1 else 20,
                  q = 0.1 * min(1, NROW(x) / NCOL(x)),
                  theta1 = 1,
                  theta2 = 0.5,
                  prox_method = c("stack", "pava"),
                  screen = TRUE,
                  screen_alg = c("strong", "previous"),
                  tol_dev_change = 1e-5,
                  tol_dev_ratio = 0.995,
                  max_variables = NROW(x),
                  solver = c("fista", "admm"),
                  max_passes = 1e6,
                  tol_abs = 1e-5,
                  tol_rel = 1e-4,
                  tol_rel_gap = 1e-5,
                  tol_infeas = 1e-3,
                  tol_rel_coef_change = 1e-3,
                  diagnostics = FALSE,
                  verbosity = 0,
                  sigma,
                  n_sigma,
                  lambda_min_ratio) {
  if (!missing(sigma)) {
    warning("`sigma` argument is deprecated; please use `alpha` instead")
    alpha <- sigma
  }

  if (!missing(n_sigma)) {
    warning(
      "`n_sigma` argument is deprecated; please use `path_length` instead"
    )
    path_length <- n_sigma
  }

  if (!missing(lambda_min_ratio)) {
    warning(
      "`lambda_min_ratio` is deprecated; please use `alpha_min_ratio` instead"
    )
    alpha_min_ratio <- lambda_min_ratio
  }

  ocall <- match.call()

  family <- match.arg(family)
  solver <- match.arg(solver)
  screen_alg <- match.arg(screen_alg)
  prox_method_choice <- switch(match.arg(prox_method), stack = 0, pava = 1)

  if (solver == "admm" && family != "gaussian") {
    stop("ADMM solver is only supported with `family = 'gaussian'`")
  }

  if (is.character(scale)) {
    scale <- match.arg(scale)
  } else if (is.logical(scale) && length(scale) == 1L) {
    scale <- ifelse(scale, "l2", "none")
  } else {
    stop("`scale` must be logical or a character")
  }

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
    tol_rel_gap >= 0,
    tol_infeas >= 0,
    tol_abs >= 0,
    tol_rel >= 0,
    is.logical(center),
    tol_rel_coef_change >= 0,
    is.numeric(tol_rel_coef_change),
    theta1 >= 0,
    theta2 >= 0,
    is.finite(theta1),
    is.finite(theta2)
  )

  fit_intercept <- intercept

  # convert sparse x to dgCMatrix class from package Matrix.
  is_sparse <- inherits(x, "sparseMatrix")

  if (NROW(y) != NROW(x)) {
    stop("the number of samples in `x` and `y` must match")
  }

  if (NROW(y) == 0) {
    stop("`y` is empty")
  }

  if (NROW(x) == 0) {
    stop("`x` is empty")
  }

  if (anyNA(y) || anyNA(x)) {
    stop("missing values are not allowed")
  }

  if (is_sparse) {
    x <- methods::as(x, "dgCMatrix")
  } else {
    x <- as.matrix(x)
  }

  if (is_sparse && center) {
    stop("centering would destroy sparsity in `x` (predictor matrix)")
  }

  res <- preprocessResponse(family, y, fit_intercept)
  y <- as.matrix(res$y)
  y_center <- res$y_center
  y_scale <- res$y_scale
  class_names <- res$class_names
  m <- n_targets <- res$n_targets
  response_names <- res$response_names
  variable_names <- colnames(x)
  max_variables <- max_variables * m

  if (is.null(variable_names)) {
    variable_names <- paste0("V", seq_len(p))
  }
  if (is.null(response_names)) {
    response_names <- paste0("y", seq_len(m))
  }

  if (is.character(alpha)) {
    alpha <- match.arg(alpha)

    if (alpha == "path") {
      alpha_type <- "auto"
      alpha <- double(path_length)
    } else if (alpha == "estimate") {
      if (family != "gaussian") {
        stop("`alpha = 'estimate'` can only be used if `family = 'gaussian'`")
      }

      alpha_type <- "estimate"
      alpha <- NULL

      if (path_length > 1) {
        warning("`path_length` ignored since `alpha = 'estimate'`")
      }
    }
  } else {
    alpha <- as.double(alpha)
    alpha_type <- "user"

    alpha <- as.double(alpha)
    path_length <- length(alpha)

    stopifnot(path_length > 0)

    if (any(alpha < 0)) {
      stop("`alpha` cannot contain negative values")
    }

    if (is.unsorted(rev(alpha))) {
      stop("`alpha` must be decreasing")
    }

    if (anyDuplicated(alpha) > 0) {
      stop("all values in `alpha` must be unique")
    }

    # do not stop path early if user requests specific alpha
    tol_dev_change <- 0
    tol_dev_ratio <- 1
    max_variables <- (NCOL(x) + intercept) * m
  }

  n_lambda <- m * p

  if (is.character(lambda)) {
    lambda_type <- match.arg(lambda)

    if (lambda_type == "bhq") {
      warning(
        "'bhq' option to argument lambda has been depracted and will",
        "will be defunct in the next release; please use 'bh' instead"
      )
    }

    lambda <- double(n_lambda)
  } else {
    lambda_type <- "user"
    lambda <- as.double(lambda)

    if (length(lambda) != m * p) {
      stop("`lambda` must be as long as there are variables")
    }

    if (is.unsorted(rev(lambda))) {
      stop("`lambda` must be non-increasing")
    }

    if (any(lambda < 0)) {
      stop("`lambda` cannot contain negative values")
    }
  }

  control <- list(
    family = family,
    fit_intercept = fit_intercept,
    is_sparse = is_sparse,
    scale = scale,
    center = center,
    path_length = path_length,
    prox_method_choice = prox_method_choice,
    n_targets = n_targets,
    screen = screen,
    screen_alg = screen_alg,
    alpha = alpha,
    alpha_type = alpha_type,
    lambda = lambda,
    lambda_type = lambda_type,
    alpha_min_ratio = alpha_min_ratio,
    q = q,
    theta1 = theta1,
    theta2 = theta2,
    y_center = y_center,
    y_scale = y_scale,
    max_passes = max_passes,
    diagnostics = diagnostics,
    verbosity = verbosity,
    max_variables = max_variables,
    solver = solver,
    tol_dev_change = tol_dev_change,
    tol_dev_ratio = tol_dev_ratio,
    tol_rel_gap = tol_rel_gap,
    tol_infeas = tol_infeas,
    tol_abs = tol_abs,
    tol_rel = tol_rel,
    tol_rel_coef_change = tol_rel_coef_change
  )

  fitSLOPE <- if (is_sparse) sparseSLOPE else denseSLOPE

  if (intercept) {
    x <- cbind(1, x)
  }

  if (alpha_type %in% c("path", "user")) {
    fit <- fitSLOPE(x, y, control)
  } else {
    # estimate the noise level, if possible
    if (is.null(alpha) && n >= p + 30) {
      alpha <- estimateNoise(x, y)
    }

    # run the solver, iteratively if necessary.
    if (is.null(alpha)) {
      # Run Algorithm 5 of Section 3.2.3. (Bogdan et al.)
      if (intercept) {
        selected <- 1
      } else {
        selected <- integer(0)
      }

      repeat {
        selected_prev <- selected

        alpha <- estimateNoise(x[, selected, drop = FALSE], y, intercept)
        control$alpha <- alpha

        fit <- fitSLOPE(x, y, control)

        selected <- which(abs(drop(fit$betas)) > 0)

        if (fit_intercept) {
          selected <- union(1, selected)
        }

        if (identical(selected, selected_prev)) {
          break
        }

        if (length(selected) + 1 >= n) {
          stop("selected >= n-1 variables; cannot estimate variance")
        }
      }
    } else {
      control$alpha <- alpha
      fit <- fitSLOPE(x, y, control)
    }
  }

  lambda <- fit$lambda
  alpha <- fit$alpha
  path_length <- length(alpha)
  active_sets <- lapply(drop(fit$active_sets), function(x) drop(x) + 1)
  beta <- fit$betas
  nonzeros <- apply(beta, c(2, 3), function(x) abs(x) > 0)
  coefficients <- beta

  if (fit_intercept) {
    nonzeros <- nonzeros[-1, , , drop = FALSE]
    dimnames(coefficients) <- list(
      c("(Intercept)", variable_names),
      response_names[1:n_targets],
      paste0("p", seq_len(path_length))
    )
  } else {
    dimnames(coefficients) <- list(
      variable_names,
      response_names[1:n_targets],
      paste0("p", seq_len(path_length))
    )
  }

  # check if maximum number of passes where reached anywhere
  passes <- fit$passes
  reached_max_passes <- passes >= max_passes

  if (any(reached_max_passes)) {
    reached_max_passes_where <- which(reached_max_passes)
    warning(
      "maximum number of passes reached at steps ",
      paste(reached_max_passes_where, collapse = ", "), "!"
    )
  }

  diagnostics <- if (diagnostics) setupDiagnostics(fit) else NULL

  slope_class <- switch(
    family,
    gaussian = "GaussianSLOPE",
    binomial = "BinomialSLOPE",
    poisson = "PoissonSLOPE",
    multinomial = "MultinomialSLOPE"
  )

  structure(
    list(
      coefficients = coefficients,
      nonzeros = nonzeros,
      lambda = lambda,
      alpha = alpha,
      class_names = class_names,
      passes = passes,
      violations = fit$violations,
      active_sets = active_sets,
      unique = drop(fit$n_unique),
      deviance_ratio = drop(fit$deviance_ratio),
      null_deviance = fit$null_deviance,
      family = family,
      diagnostics = diagnostics,
      call = ocall
    ),
    class = c(slope_class, "SLOPE")
  )
}
