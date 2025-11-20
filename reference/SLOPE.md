# Sorted L-One Penalized Estimation

Fit a generalized linear model regularized with the sorted L1 norm,
which applies a non-increasing regularization sequence to the
coefficient vector (\\\beta\\) after having sorted it in decreasing
order according to its absolute values.

## Usage

``` r
SLOPE(
  x,
  y,
  family = c("gaussian", "binomial", "multinomial", "poisson"),
  intercept = TRUE,
  center = c("mean", "min", "none"),
  scale = c("sd", "l1", "l2", "max_abs", "none"),
  alpha = c("path", "estimate"),
  lambda = c("bh", "gaussian", "oscar", "lasso"),
  alpha_min_ratio = if (NROW(x) < NCOL(x)) 0.01 else 1e-04,
  path_length = 100,
  q = 0.1,
  theta1 = 1,
  theta2 = 0.5,
  tol_dev_change = 1e-05,
  tol_dev_ratio = 0.999,
  max_variables = NROW(x) + 1,
  solver = c("auto", "hybrid", "pgd", "fista", "admm"),
  max_passes = 1e+06,
  tol = 1e-04,
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
)
```

## Arguments

- x:

  the design matrix, which can be either a dense matrix of the standard
  *matrix* class, or a sparse matrix inheriting from
  [Matrix::sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).
  Data frames will be converted to matrices internally.

- y:

  the response, which for `family = "gaussian"` must be numeric; for
  `family = "binomial"` or `family = "multinomial"`, it can be a factor.

- family:

  model family (objective); see **Families** for details.

- intercept:

  whether to fit an intercept

- center:

  whether to center predictors or not by their mean. Defaults to `TRUE`
  if `x` is dense and `FALSE` otherwise.

- scale:

  type of scaling to apply to predictors.

  - `"l1"` scales predictors to have L1 norms of one.

  - `"l2"` scales predictors to have L2 norms of one.#'

  - `"sd"` scales predictors to have a population standard deviation
    one.

  - `"none"` applies no scaling.

- alpha:

  scale for regularization path: either a decreasing numeric vector
  (possibly of length 1) or a character vector; in the latter case, the
  choices are:

  - `"path"`, which computes a regularization sequence where the first
    value corresponds to the intercept-only (null) model and the last to
    the almost-saturated model, and

  - `"estimate"`, which estimates a *single* `alpha` using Algorithm 5
    in Bogdan et al. (2015).

  When a value is manually entered for `alpha`, it will be scaled based
  on the type of standardization that is applied to `x`. For
  `scale = "l2"`, `alpha` will be scaled by \\\sqrt n\\. For
  `scale = "sd"` or `"none"`, alpha will be scaled by \\n\\, and for
  `scale = "l1"` no scaling is applied. Note, however, that the `alpha`
  that is returned in the resulting value is the **unstandardized**
  alpha.

- lambda:

  either a character vector indicating the method used to construct the
  lambda path or a numeric non-decreasing vector with length equal to
  the number of coefficients in the model; see section **Regularization
  sequences** for details.

- alpha_min_ratio:

  smallest value for `lambda` as a fraction of `lambda_max`; used in the
  selection of `alpha` when `alpha = "path"`.

- path_length:

  length of regularization path; note that the path returned may still
  be shorter due to the early termination criteria given by
  `tol_dev_change`, `tol_dev_ratio`, and `max_variables`.

- q:

  parameter controlling the shape of the lambda sequence, with usage
  varying depending on the type of path used and has no effect is a
  custom `lambda` sequence is used. Must be greater than `1e-6` and
  smaller than 1.

- theta1:

  parameter controlling the shape of the lambda sequence when
  `lambda == "OSCAR"`. This parameter basically sets the intercept for
  the lambda sequence and is equivalent to \\\lambda_1\\ in the original
  OSCAR formulation.

- theta2:

  parameter controlling the shape of the lambda sequence when
  `lambda == "OSCAR"`. This parameter basically sets the slope for the
  lambda sequence and is equivalent to \\\lambda_2\\ in the original
  OSCAR formulation.

- tol_dev_change:

  the regularization path is stopped if the fractional change in
  deviance falls below this value; note that this is automatically set
  to 0 if a alpha is manually entered

- tol_dev_ratio:

  the regularization path is stopped if the deviance ratio \\1 -
  \mathrm{deviance}/\mathrm{(null-deviance)} \\ is above this threshold

- max_variables:

  criterion for stopping the path in terms of the maximum number of
  unique, nonzero coefficients in absolute value in model. For the
  multinomial family, this value will be multiplied internally with the
  number of levels of the response minus one.

- solver:

  type of solver use, either `"auto"`, `"hybrid"`, `"pgd"`, or
  `"fista"`; `"auto"` means that the solver is automatically selected,
  which currently means that `"hybrid"` is used for all objectives
  except multinomial ones, in which case FISTA (`"fista"`) is used.

- max_passes:

  maximum number of passes (outer iterations) for solver

- tol:

  stopping criterion for the solvers in terms of the relative duality
  gap

- threads:

  number of threads to use in the solver; if `NULL`, half of the
  available (logical) threads will be used

- diagnostics:

  whether to save diagnostics from the solver (timings and other values
  depending on type of solver)

- patterns:

  whether to return the SLOPE pattern (cluster, ordering, and sign
  information) as a list of sparse matrices, one for each step on the
  path.

- gamma:

  relaxation mixing parameter, between 0 and 1. Has no effect if set
  to 0. If larger than 0, the solver will mix the coefficients from the
  ordinary SLOPE solutions with the coefficients from the relaxed
  solutions (fitting OLS on the SLOPE pattern).

- cd_type:

  Type of coordinate descent to use, either `"cyclical"` or
  `"permuted"`. The former means that the cluster are cycled through in
  descending order of their coefficients' magnitudes, while the latter
  means that the clusters are permuted randomly for each pass.

- tol_abs:

  DEPRECATED

- tol_rel:

  relative DEPRECATED

- tol_rel_gap:

  DEPRECATED

- tol_infeas:

  DEPRECATED

- tol_rel_coef_change:

  DEPRECATED

- prox_method:

  DEPRECATED

- screen:

  DEPRECATED

- verbosity:

  DEPRECATED

- screen_alg:

  DEPRECATED

## Value

An object of class `"SLOPE"` with the following slots:

- coefficients:

  a list of the coefficients from the model fit, not including the
  intercept. The coefficients are stored as sparse matrices.

- nonzeros:

  a three-dimensional logical array indicating whether a coefficient was
  zero or not

- lambda:

  the lambda vector that when multiplied by a value in `alpha` gives the
  penalty vector at that point along the regularization path

- alpha:

  vector giving the (unstandardized) scaling of the lambda sequence

- class_names:

  a character vector giving the names of the classes for binomial and
  multinomial families

- passes:

  the number of passes the solver took at each step on the path

- deviance_ratio:

  the deviance ratio (as a fraction of 1)

- null_deviance:

  the deviance of the null (intercept-only) model

- family:

  the name of the family used in the model fit

- diagnostics:

  a `data.frame` of objective values for the primal and dual problems,
  as well as a measure of the infeasibility, time, and iteration; only
  available if `diagnostics = TRUE` in the call to `SLOPE()`.

- n_observations:

  the number of observations in the training data

- n_predictors:

  the number of predictors in the training data

- call:

  the call used for fitting the model

## Details

`SLOPE()` solves the convex minimization problem \$\$ f(\beta) + \alpha
\sum\_{i=j}^p \lambda_j \|\beta\|\_{(j)}, \$\$ where \\f(\beta)\\ is a
smooth and convex function and the second part is the sorted L1-norm. In
ordinary least-squares regression, \\f(\beta)\\ is simply the squared
norm of the least-squares residuals. See section **Families** for
specifics regarding the various types of \\f(\beta)\\ (model families)
that are allowed in `SLOPE()`.

By default, `SLOPE()` fits a path of models, each corresponding to a
separate regularization sequence, starting from the null
(intercept-only) model to an almost completely unregularized model.
These regularization sequences are parameterized using \\\lambda\\ and
\\\alpha\\, with only \\\alpha\\ varying along the path. The length of
the path can be manually, but will terminate prematurely depending on
arguments `tol_dev_change`, `tol_dev_ratio`, and `max_variables`. This
means that unless these arguments are modified, the path is not
guaranteed to be of length `path_length`.

## Families

**Gaussian**

The Gaussian model (Ordinary Least Squares) minimizes the following
objective: \$\$ \frac{1}{2} \Vert y - X\beta\Vert_2^2 \$\$

**Binomial**

The binomial model (logistic regression) has the following objective:
\$\$ \sum\_{i=1}^n \log\left(1+ \exp\left( - y_i \left(x_i^T\beta +
\beta_0 \right) \right) \right) \$\$ with \\y \in \\-1, 1\\\\.

**Poisson**

In poisson regression, we use the following objective:

\$\$ -\sum\_{i=1}^n \left(y_i\left( x_i^T\beta + \beta_0\right) -
\exp\left(x_i^T\beta + \beta_0 \right)\right) \$\$

**Multinomial**

In multinomial regression, we minimize the full-rank objective \$\$
-\sum\_{i=1}^n\left( \sum\_{k=1}^{m-1} y\_{ik}(x_i^T\beta_k +
\beta\_{0,k}) - \log\sum\_{k=1}^{m-1} \exp\big(x_i^T\beta_k +
\beta\_{0,k}\big) \right) \$\$ with \\y\_{ik}\\ being the element in a
\\n\\ by \\(m-1)\\ matrix, where \\m\\ is the number of classes in the
response.

## Regularization Sequences

There are multiple ways of specifying the `lambda` sequence in
`SLOPE()`. It is, first of all, possible to select the sequence manually
by using a non-increasing numeric vector, possibly of length one, as
argument instead of a character. The greater the differences are between
consecutive values along the sequence, the more clustering behavior will
the model exhibit. Note, also, that the scale of the \\\lambda\\ vector
makes no difference if `alpha = NULL`, since `alpha` will be selected
automatically to ensure that the model is completely sparse at the
beginning and almost unregularized at the end. If, however, both `alpha`
and `lambda` are manually specified, then the scales of both do matter,
so make sure to choose them wisely.

Instead of choosing the sequence manually, one of the following
automatically generated sequences may be chosen.

**BH (Benjamini–Hochberg)**

If `lambda = "bh"`, the sequence used is that referred to as
\\\lambda^{(\mathrm{BH})}\\ by Bogdan et al, which sets \\\lambda\\
according to \$\$ \lambda_i = \Phi^{-1}(1 - iq/(2p)), \$\$ for
\\i=1,\dots,p\\, where \\\Phi^{-1}\\ is the quantile function for the
standard normal distribution and \\q\\ is a parameter that can be set by
the user in the call to `SLOPE()`.

**Gaussian**

This penalty sequence is related to BH, such that \$\$ \lambda_i =
\lambda^{(\mathrm{BH})}\_i \sqrt{1 + w(i-1)\cdot
\mathrm{cumsum}(\lambda^2)\_i}, \$\$ for \\i=1,\dots,p\\, where \\w(k) =
1/(n-k-1)\\. We let \\\lambda_1 = \lambda^{(\mathrm{BH})}\_1\\ and
adjust the sequence to make sure that it's non-increasing. Note that if
\\p\\ is large relative to \\n\\, this option will result in a constant
sequence, which is usually not what you would want.

**OSCAR**

This sequence comes from Bondell and Reich and is a linear
non-increasing sequence, such that \$\$ \lambda_i = \theta_1 + (p -
i)\theta_2. \$\$ for \\i = 1,\dots,p\\. We use the parametrization from
Zhong and Kwok (2021) but use \\\theta_1\\ and \\\theta_2\\ instead of
\\\lambda_1\\ and \\\lambda_2\\ to avoid confusion and abuse of
notation.

**lasso**

SLOPE is exactly equivalent to the lasso when the sequence of
regularization weights is constant, i.e. \$\$ \lambda_i = 1 \$\$ for \\i
= 1,\dots,p\\. Here, again, we stress that the fact that all \\\lambda\\
are equal to one does not matter as long as `alpha == NULL` since we
scale the vector automatically. Note that this option is only here for
academic interest and to highlight the fact that SLOPE is a
generalization of the lasso. There are more efficient packages, such as
**glmnet** and **biglasso**, for fitting the lasso.

## Solvers

There are currently three solvers available for SLOPE: Hybrid (Beck and
Teboulle 2009), proximal gradient descent (PGD), and FISTA (Beck and
Teboulle, 2009). The hybrid method is the preferred and generally
fastest method and is therefore the default for the Gaussian and
binomial families, but not currently available for multinomial and
disabled for Poisson due to convergence issues.

## References

Bogdan, M., van den Berg, E., Sabatti, C., Su, W., & Candès, E. J.
(2015). SLOPE – adaptive variable selection via convex optimization. The
Annals of Applied Statistics, 9(3), 1103–1140.

Larsson, J., Klopfenstein, Q., Massias, M., & Wallin, J. (2023).
Coordinate descent for SLOPE. In F. Ruiz, J. Dy, & J.-W. van de Meent
(Eds.), Proceedings of the 26th international conference on artificial
intelligence and statistics (Vol. 206, pp. 4802–4821). PMLR.
https://proceedings.mlr.press/v206/larsson23a.html

Bondell, H. D., & Reich, B. J. (2008). Simultaneous Regression
Shrinkage, Variable Selection, and Supervised Clustering of Predictors
with OSCAR. Biometrics, 64(1), 115–123. JSTOR.

Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2010).
Distributed Optimization and Statistical Learning via the Alternating
Direction Method of Multipliers. Foundations and Trends® in Machine
Learning, 3(1), 1–122.

Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding
Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences,
2(1), 183–202.

## See also

[`plot.SLOPE()`](https://jolars.github.io/SLOPE/reference/plot.SLOPE.md),
[`plotDiagnostics()`](https://jolars.github.io/SLOPE/reference/plotDiagnostics.md),
[`score()`](https://jolars.github.io/SLOPE/reference/score.md),
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md),
[`trainSLOPE()`](https://jolars.github.io/SLOPE/reference/trainSLOPE.md),
[`coef.SLOPE()`](https://jolars.github.io/SLOPE/reference/coef.SLOPE.md),
[`print.SLOPE()`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md),
[`print.SLOPE()`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md),
[`deviance.SLOPE()`](https://jolars.github.io/SLOPE/reference/deviance.SLOPE.md),
[`sortedL1Prox()`](https://jolars.github.io/SLOPE/reference/sortedL1Prox.md)

## Examples

``` r
# Gaussian response, default lambda sequence
fit <- SLOPE(bodyfat$x, bodyfat$y)

# Multinomial response, custom alpha and lambda
m <- length(unique(wine$y)) - 1
p <- ncol(wine$x)

alpha <- 0.005
lambda <- exp(seq(log(2), log(1.8), length.out = p * m))

fit <- SLOPE(
  wine$x,
  wine$y,
  family = "multinomial",
  lambda = lambda,
  alpha = alpha
)
```
