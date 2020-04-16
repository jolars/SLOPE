# SLOPE 0.2.1

## Minor changes

* A few examples in `deviance()` and `SLOPE()` that were taking
  too long to execute have been removed or modified.

# SLOPE 0.2.0

This version of SLOPE represents a major change to the package. We have
merged functionality from the owl package into this package, which
means there are several changes to the API, including deprecated functions.

## Major changes

* `SLOPE_solver()`, `SLOPE_solver_matlab()`, `prox_sorted_L1()`,
  and `create_lambda()`
  have been deprecated (and will be defunct in the
  next version of SLOPE)
* arguments `X`, `fdr`, and `normalize` have been deprecated
  in `SLOPE()` and replaced by `x`, `q`, `scale` and `center`, respectively
* options `"default"` and `"matlab"` to argument
  `solver` in `SLOPE()` have been deprecated and replaced with `"fista"`
  and `"admm"`, which uses the accelerated proximal gradient method
  FISTA and alternating direction of multipliers method (ADMM)
  respectively
* ADMM has been implemented as a solver for `family = "gaussian"`
* binomial, poisson, and multinomial families are now supported (using
  `family` argument in `SLOPE()`)
* input to `lambda` is now scaled (divided by) the number of observations (rows)
  in `x`
* predictor screening rules have been implemented and are activated by
  calling `screen = TRUE` in `SLOPE()`. The type of algorithm can also
  be set via `screen_alg`.
* `SLOPE()` now returns an object of class `"SLOPE"` (and an additional
  class depending on input to `family` in `SLOPE()`
* `SLOPE` objects gain `coef()` and `plot()` methods.
* `SLOPE` now uses screening rules to speed up model fitting in the
  high-dimensional regime
* most of the code is now written in C++ using the **Rcpp** and **RcppArmadillo**
  packages
* a new function `trainSLOPE()` trains SLOPE with repeated k-folds 
  cross-validation
* a new function `caretSLOPE()` enables model-tuning using the
  **caret** package
* `SLOPE()` now fits an entire path of regularization sequences by default
* the `normalize` option to `SLOPE()` has been replaced by `scale` and
  `center`, which allows granular options for standardization
* sparse matrices (from the **Matrix** package) can now be used as
  input
* there are now five datasets included in the package
* the introductory vignette has been replaced

## Minor changes

* a new function `deviance()` returns the deviance from the fit
* a new function `score()` can be used to assess model performance against
  new data
* a new function `plotDiagnostics()` has been included to visualize
  data from the solver (if `diagnostics = TRUE` in the call to `SLOPE()`)
* OSCAR-type penalty sequences can be used by setting `lambda = "oscar"
  in the call to `SLOPE()`
* the test suite for the package has been greatly extended

