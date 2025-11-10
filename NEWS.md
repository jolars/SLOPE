# SLOPE 1.2.0

## Major Changes

- The default argument for `threads` in `SLOPE()` has been changed from `NULL`
  (half of available cores) to `1` (no multithreading). This is to avoid
  excessive CPU usage on systems with many cores. Users can still set
  `threads = NULL` to get the previous behavior, or set it to any positive
  integer to control the number of threads used.

## Minor Changes

- The C++ routine can now be interrupted through R's usual interrupt mechanism
  (e.g., Ctrl+C in RStudio or R terminal).
- The citation information has been updated. Use `citation("SLOPE")` to see the
  correct way to cite the package.
- The y axis label for coefficients has changed from `expression(hat(beta))` to
  `"Coefficients"` to avoid issues with rendering in some environments, where
  the hat symbol would be cropped.

## Bug Fixes

- In the return object from `SLOPE()`, the `coefficients_scaled` field
  incorrectly contained unscaled coefficients, which also affected
  `coef.SLOPE()`. This has now been fixed.

# SLOPE 1.1.0

## New Features

- A new data set, `glioma`, has been added to the package. It contains gene
  expression measurements for patients with glioma along with healthy controls.

# SLOPE 1.0.1

## Bug Fixes

- Fixed a test error on the M1mac platform.

# SLOPE 1.0.0

This update of SLOPE brings an entirely different C++ implementation of the
underlying package based on the C++ library
[libslope](https://github.com/jolars/libslope). It comes with several large and
breaking changes with respect to the previous version of the package.

We realized that this may throw off some users, and hope that you will be
patient with dealing with the large number of breaking changes.

## Breaking Changes

- The `caretSLOPE()` function that was deprecated has now been removed from the
  package.
- Fields `unique`, `violations`, and `active_sets` are no longer stored in the
  `SLOPE` object. These fields were typically only used for debugging purposes.
- The `prox_method` and `method` arguments in `SLOPE()` and `sortedL1Prox()`,
  respectively, have been removed. The proximal operator is now always computed
  using the fast stack-based algorithm. There was never any reason to use the
  slower PAVA algorithm.
- The ADMM solver has been removed from the package. Calling `SLOPE()` with
  `solver = "admm"` will now throws a warning and the value will be
  automatically set to `"auto"`.
- `alpha` is now scaled by `n` (the number of observations) and differences with
  respect to the type of scaling are no longer taken into account.
- The object `coefficients` from `SLOPE()` is now a list of sparse matrices
  (rather than a three-dimensional array as before). Now it contains only the
  coefficients and not the intercepts. The intercepts are instead stored in
  `intercepts` in the returned object and are always present even if
  `intercept = FALSE`.
- The behavior of `coef.SLOPE()` has changed somewhat, and if
  `simplify = FALSE`, then the returned object is now instead a list of sparse
  matrices (rather than a three-dimensional array as before).
- The default value of `q` in `SLOPE()` has changed from
  `0.1 * min(1, NROW(x) / NCOL(x))` to `0.1`.
- Arguments `sigma`, `n_sigma`, and `lambda_min_ratio` in `SLOPE()` that were
  previously deprecated have been removed.
- `SLOPE()` now internally solves the problem normalized by scaling with the
  number of observations, which means that values returned in `deviance` and
  `prmals` and `duals` if `diagnostics = TRUE` are now scaled by `n`.
- `path_length` in `SLOPE()` now defaults to 100 (previously 20).
- `tol_dev_ratio` in `SLOPE()` now defaults to `0.999` (previously `0.995`).
- Plots from `plot.SLOPE()` now use base R graphics rather than ggplot2. This
  means that the plots are more difficult to customize but plot much more faster
  when there are many variables and significantly reduces the dependency load of
  the package. For plots of trained SLOPE objects, which used to be faceted on
  the `q` parameter, the user now needs to use the standard base R graphics API
  to facet plots via `par(mfrow = c(1, 2))` or similar.

## Deprecated Functionality

- Arguments `tol_rel_gap`, `tol_infeas`, `tol_abs`, `tol_rel`, `tol_rel_coef` in
  `SLOPE()` are now deprecated. The solvers now all rely on the same tolerance
  criterion, which is set by `tol` and uses the duality gap normalized by the
  current primal value.
- Arguments `screen` and `screen_alg` are now deprecated and have no effect.
  Feature screening is always used. These arguments were only used for
  debugging.
- The argument `verbosity` in `SLOPE()` is now defunct and has no effect.
- The argument `prox_method` in `SLOPE()` and `sortedL1Prox()` is now defunct
  and has no effect.

## New Features

- Centering `x` in `SLOPE()` is now allowed again, even when the matrix is
  sparse.
- Out-of-memory matrices are now allowed through the `bigmemory` package. Only
  support for dense matrices is available at the moment.
- Centers and scales can now be specified manually by providing vectors to
  `center` and `scale` in `SLOPE()`.
- A new solver based on a hybrid method of proximal gradient descent and
  coordinate descent is available and used by default by the Gaussian and
  binomial families. Use it by specifying `solver = "hybrid"`.
- Solver can now be set to `"auto"`, in which case the package automatically
  chooses a solver.
- The returned duality gaps when `diagnostics = TRUE` are now _true_ duality
  gaps, computed by guaranteeing that the dual variable is feasible (which was
  not the case previously).
- `scale` in `SLOPE()` gains a new option `"max_abs"` which scales the columns
  of `x` by their maximum absolute value.
- When `alpha = "estimate"`, there is a now an iteration limit in case the
  algorithm does not converge to one set of features. Thanks @RomanParzer.
- `plot.SLOPE()` gains a new argument `magnitudes`, which causes the plot to
  only show the magnitudes of the coefficients (which helps if you want to
  visualize cluster structure).
- `plot.SLOPE()` gains a new argument `add_labels`, which add numbers for the
  coefficients to the plot. Set to `FALSE` by default.
- Relaxed SLOPE models can now be fit by specifying `gamma` in `SLOPE()`.
- `plot.trainedSLOPE()` gains a new argument `index`, to select which of the
  hyperparameter combinations to plot for.
- There's a new function `plotClusters()`, which allows plotting the cluster
  structure in SLOPE. Thanks, @KrystynaGrzesiak!
- `SLOPE()` gains a new argument `cd_type`, to control the type of coordinate
  descent used for the hybrid solver, with options `"cyclical"` and
  `"permuted"`.

## Bug Fixes

- Return correct model when training for AUC in `trainSLOPE()`.

## Performance Improvements

The new hybrid algorithm that's implemented in libslope and now used in the
package constitutes a major upgrade in terms of performance.

- The solver is now much more memory-efficient and can avoid copies of the
  design matrix entirely by normalizing the columns just-in-time. This is the
  standard behavior. Future versions of the package will allow the user to
  specify whether to copy (and modify) the design matrix or not.

## Dependencies

We have made an effort to reduce the footprint of the package and reduce the
number of dependencies.

- The package now relies on Eigen (through RcppEigen) rather than Armadillo,
  which means that there is no longer any reliance on BLAS and LAPACK libraries.
- The dependency on `ggplot2` is removed.
- The `vdiffr`, `tidyr`, `dplyr`, `bench`, `scales`, and `glmnet` packages in
  the `Suggests` field that were used for testing are now removed.

# SLOPE 0.5.2

## Bug Fixes

- Fixed bug when computing regularization weights for type `"gaussian"` when the
  number of observations is less than the number of variables.

# SLOPE 0.5.1

## Minor Changes

- Website updated to bootstrap 5-based pkgdown theme.
- Updated e-mail of maintainer.
- Dependencies on checkmate and mice were dropped.
- Update sparse matrix coercion to avoid deprecated functionality in the Matrix
  package.

# SLOPE 0.5.0

## Major changes

- `plot.SLOPE()`, `plot.trainSLOPE()` and `plotDiagnostics()` have been
  reimplemented in ggplot2.

## Deprecated Functions

- `caretSLOPE()` has been deprecated and will be made defunct in version 0.6.0.

# SLOPE 0.4.1

## Bug Fixes

- The C++ standard library _memory_ was added to a source file to fix
  compilation errors on some systems.

# SLOPE 0.4.0

## New Functions

- `sortedL1Prox()` is a new function that computes the proximal operator for the
  sorted L1 norm (the penalty term in SLOPE).
- `regularizationWeights()` is a new function that returns the penalty weights
  (lambda sequence) for SLOPE or OSCAR.

## Major changes

- The parametrization for OSCAR models have been corrected and changed. As a
  result, `SLOPE()` gains two arguments: `theta1` and `theta2` to control the
  behavior using the parametrization from L. W. Zhong and J. T. Kwok, “Efficient
  sparse modeling with automatic feature grouping,” IEEE Transactions on Neural
  Networks and Learning Systems, vol. 23, no. 9, pp. 1436–1447, Sep. 2012, doi:
  10.1109/TNNLS.2012.2200262. `q` is no longer used with OSCAR models. Thanks,
  Nuno Eusebio.
- `SLOPE()` has gained a new argument, `prox_method`, which allows the user to
  select prox algorithm to use. There is no an additional algorithm in the
  package, based on the PAVA algorithm used in isotonic regression, that can be
  used. Note that this addition is mostly of academic interest and does not need
  to be changed by the user.

## Minor Changes

- The `q` parameter is no longer allowed to be smaller than `1e-6` to avoid
  constructions of regularization paths with infinite `lambda` values.
- The `lambda` argument in `SLOPE()` now also allowed the input `"lasso"` to
  obtain the standard lasso.
- The performance of `trainSLOPE()`

## Vignettes

- A new vignette has been added to compare algorithms for the proximal operator.

## Bug Fixes

- For very small numbers of observations (10 or so), the regularization weights
  for `lambda = "gaussian"` were incorrectly computed, increasing and then
  decreasing. This is now fixed and regularization weights in this case are now
  always non-increasing.
- Misclassification error was previously computed incorrectly in `trainSLOPE()`
  for multinomial models (thanks @jakubkala and @KrystynaGrzesiak)
- Performance of `trainSLOPE()` was previously hampered by erroneous refitting
  of the models, which has been fixed now (thanks @jakubkala and
  @KrystynaGrzesiak)

## Deprecated and Defunct

- `yvar` argument in `plotDiagnostics()` that was previously deprecated is now
  defunct.
- Using `missclass` for the `measure` argument in `trainSLOPE()` has been
  deprecated in favor of `misclass`.

# SLOPE 0.3.3

## Bug fixes

- Fixed first coefficient missing from plot if no intercept was used in the call
  to `SLOPE()`.
- Fixed incorrect results when `intercept = FALSE` and `family = "gaussian"`
  (#13, thanks, Patrick Tardivel).

# SLOPE 0.3.2

## Minor changes

- Added `tol_rel_coef_change` argument to `SLOPE()` as a convergence criterion
  for the FISTA solver that sets a tolerance for the relative change in
  coefficients across iterations.

## Bug fixes

- Fixed premature stopping of the solver for the first step of the
  regularization path (the null model).
- Actually fix UBSAN/ASAN sanitizer warnings by modifying code for FISTA solver.

# SLOPE 0.3.1

## Bug fixes

- Fixed package build breaking on solaris because of missing STL namespace
  specifier for `std::sqrt()` in `src/SLOPE.cpp`.
- Fixed erroneous scaling of absolute tolerance in stopping criteria for the
  ADMM solvers. Thanks, @straw-boy.
- Fixed sanitizer warning from CRAN checks.

# SLOPE 0.3.0

## Major changes

- Scaling of `alpha` (previously `sigma`) is now invariant to the number of
  observations, which is achieved by scaling the penalty part of the objective
  by the square root of the number of observations if `scale = "l2"` and the
  number of observations if `scale = "sd"` or `"none"`. No scaling is applied
  when `scale = "l1"`.
- The `sigma` argument is deprecated in favor of `alpha` in `SLOPE()`,
  `coef.SLOPE()`, and `predict.SLOPE()`.
- The `n_sigma` argument is deprecated in favor of `path_length` in `SLOPE()`
- The `lambda_min_ratio` argument is deprecated in favor of `alpha_min_ratio` in
  `SLOPE()`
- The default for argument `lambda` in `SLOPE()` has changed from `"gaussian"`
  to `"bh"`.
- Functions and arguments deprecated in 0.2.0 are now defunct and have been
  removed from the package.
- `scale = "sd"` now scales with the population standard deviation rather than
  the sample standard deviation, i.e. the scaling factor now used is the number
  of observations (and not the number of observations minus one as before).

## Minor changes

- Default `path_length` has changed from 100 to 20.
- `plot.SLOPE()` has gained an argument `x_variable` that controls what is
  plotted on the x axis.
- A warning is now thrown if the maximum number of passes was reached anywhere
  along the path (and prints where as well).
- If the `max_variables` criterion is hit, the solution path returned will now
  include also the last solution (which was not the case before). Thanks,
  @straw-boy.

## Bug fixes

- Plotting models that are completely sparse no longer throws an error.
- `rho` instead of `1` is now used in the factorization part for the ADMM
  solver.

# SLOPE 0.2.1

## Minor changes

- A few examples in `deviance()` and `SLOPE()` that were taking too long to
  execute have been removed or modified.

# SLOPE 0.2.0

This version of SLOPE represents a major change to the package. We have merged
functionality from the owl package into this package, which means there are
several changes to the API, including deprecated functions.

## Major changes

- `SLOPE_solver()`, `SLOPE_solver_matlab()`, `prox_sorted_L1()`, and
  `create_lambda()` have been deprecated (and will be defunct in the next
  version of SLOPE)
- arguments `X`, `fdr`, and `normalize` have been deprecated in `SLOPE()` and
  replaced by `x`, `q`, `scale` and `center`, respectively
- options `"default"` and `"matlab"` to argument `solver` in `SLOPE()` have been
  deprecated and replaced with `"fista"` and `"admm"`, which uses the
  accelerated proximal gradient method FISTA and alternating direction of
  multipliers method (ADMM) respectively
- ADMM has been implemented as a solver for `family = "gaussian"`
- binomial, poisson, and multinomial families are now supported (using `family`
  argument in `SLOPE()`)
- input to `lambda` is now scaled (divided by) the number of observations (rows)
  in `x`
- predictor screening rules have been implemented and are activated by calling
  `screen = TRUE` in `SLOPE()`. The type of algorithm can also be set via
  `screen_alg`.
- `SLOPE()` now returns an object of class `"SLOPE"` (and an additional class
  depending on input to `family` in `SLOPE()`
- `SLOPE` objects gain `coef()` and `plot()` methods.
- `SLOPE` now uses screening rules to speed up model fitting in the
  high-dimensional regime
- most of the code is now written in C++ using the **Rcpp** and
  **RcppArmadillo** packages
- a new function `trainSLOPE()` trains SLOPE with repeated k-folds
  cross-validation
- a new function `caretSLOPE()` enables model-tuning using the **caret** package
- `SLOPE()` now fits an entire path of regularization sequences by default
- the `normalize` option to `SLOPE()` has been replaced by `scale` and `center`,
  which allows granular options for standardization
- sparse matrices (from the **Matrix** package) can now be used as input
- there are now five datasets included in the package
- the introductory vignette has been replaced

## Minor changes

- a new function `deviance()` returns the deviance from the fit
- a new function `score()` can be used to assess model performance against new
  data
- a new function `plotDiagnostics()` has been included to visualize data from
  the solver (if `diagnostics = TRUE` in the call to `SLOPE()`)
- OSCAR-type penalty sequences can be used by setting
  `lambda = "oscar" in the call to`SLOPE()`
- the test suite for the package has been greatly extended
