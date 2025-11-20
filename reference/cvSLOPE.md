# Tune SLOPE with Cross-Validation

This function trains a model fit by
[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md) by tuning
its parameters through cross-validation.

## Usage

``` r
cvSLOPE(
  x,
  y,
  q = 0.2,
  gamma = 0,
  n_folds = 10,
  n_repeats = 1,
  measure = c("mse", "mae", "deviance", "misclass", "auc"),
  refit = TRUE,
  ...
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

- q:

  a vector of quantiles for the `q` parameter in SLOPE

- gamma:

  relaxation parameter for SLOPE. Default is `0.0`, which implies to
  relaxation of the penalty.

- n_folds:

  number of folds (cross-validation)

- n_repeats:

  number of folds (cross-validation)

- measure:

  DEPRECATED

- refit:

  logical; if `TRUE`, refits the model on the full dataset using the
  optimal parameters. Default is `TRUE`.

- ...:

  other arguments to pass on to
  [`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md)

## Value

An object of class `"TrainedSLOPE"`, with the following slots:

- summary:

  a summary of the results with means, standard errors, and 0.95
  confidence levels

- data:

  the raw data from the model training

- optima:

  a `data.frame` of the best (mean) values for the different metrics and
  their corresponding parameter values

- measure:

  a `data.frame` listing the used metric and its label

- model:

  the model fit to the entire dataset using optimal parameters (only
  present if `refit = TRUE`)

- call:

  the call

## See also

Other model-tuning:
[`plot.TrainedSLOPE()`](https://jolars.github.io/SLOPE/reference/plot.TrainedSLOPE.md),
[`refit()`](https://jolars.github.io/SLOPE/reference/refit.md),
[`summary.TrainedSLOPE()`](https://jolars.github.io/SLOPE/reference/summary.TrainedSLOPE.md),
[`trainSLOPE()`](https://jolars.github.io/SLOPE/reference/trainSLOPE.md)

## Examples

``` r
# 8-fold cross-validation
tune <- cvSLOPE(
  subset(mtcars, select = c("mpg", "drat", "wt")),
  mtcars$hp,
  q = c(0.1, 0.2),
  n_folds = 8,
  n_repeats = 2,
  measure = "mse"
)

# Access the refitted model
tune$model
#> 
#> Call:
#> (function (x, y, family = c("gaussian", "binomial", "multinomial", 
#>     "poisson"), intercept = TRUE, center = c("mean", "min", "none"), 
#>     scale = c("sd", "l1", "l2", "max_abs", "none"), alpha = c("path", 
#>         "estimate"), lambda = c("bh", "gaussian", "oscar", "lasso"), 
#>     alpha_min_ratio = if (NROW(x) < NCOL(x)) 0.01 else 1e-04, 
#>     path_length = 100, q = 0.1, theta1 = 1, theta2 = 0.5, tol_dev_change = 1e-05, 
#>     tol_dev_ratio = 0.999, max_variables = NROW(x) + 1, solver = c("auto", 
#>         "hybrid", "pgd", "fista", "admm"), max_passes = 1e+06, 
#>     tol = 1e-04, threads = 1, diagnostics = FALSE, patterns = FALSE, 
#>     gamma = 1, cd_type = c("permuted", "cyclical"), tol_abs, 
#>     tol_rel, tol_rel_gap, tol_infeas, tol_rel_coef_change, prox_method, 
#>     screen, verbosity, screen_alg) 
#> {
#>     if (!missing(tol_abs)) {
#>         warning("`tol_abs` argument is deprecated and has no effect")
#>     }
#>     if (!missing(tol_rel)) {
#>         warning("`tol_rel` argument is deprecated and has no effect")
#>     }
#>     if (!missing(tol_rel_gap)) {
#>         warning("`tol_rel_gap` argument is deprecated and has no effect")
#>     }
#>     if (!missing(tol_infeas)) {
#>         warning("`tol_infeas` argument is deprecated and has no effect")
#>     }
#>     if (!missing(tol_rel_coef_change)) {
#>         warning("`tol_rel_coef_change` argument is deprecated and has no effect")
#>     }
#>     if (!missing(screen)) {
#>         warning("`screen` argument is deprecated and has no effect")
#>     }
#>     if (!missing(screen_alg)) {
#>         warning("`screen_alg` argument is deprecated and has no effect")
#>     }
#>     if (!missing(verbosity)) {
#>         warning("`verbosity` argument is deprecated and has no effect")
#>     }
#>     if (!missing(prox_method)) {
#>         warning("The 'method' argument is deprecated and ", "has no effect. It will be removed in a future version.")
#>     }
#>     ocall <- match.call()
#>     control <- processSlopeArgs(x, y, family, intercept, center, 
#>         scale, alpha, lambda, alpha_min_ratio, path_length, q, 
#>         theta1, theta2, tol_dev_change, tol_dev_ratio, max_variables, 
#>         solver, max_passes, tol, diagnostics, patterns, threads, 
#>         gamma, cd_type)
#>     x <- control$x
#>     y <- control$y
#>     alpha <- control$alpha
#>     is_sparse <- inherits(x, "sparseMatrix")
#>     is_big_matrix <- inherits(x, "big.matrix")
#>     fitSLOPE <- if (is_sparse) {
#>         sparseSLOPE
#>     }
#>     else if (is_big_matrix) {
#>         bigSLOPE
#>     }
#>     else {
#>         denseSLOPE
#>     }
#>     if (is_big_matrix) {
#>         x <- x@address
#>     }
#>     fit <- fitSLOPE(x, y, control)
#>     lambda <- fit$lambda
#>     alpha <- fit$alpha
#>     path_length <- length(alpha)
#>     intercepts <- fit$intercepts
#>     intercepts_scaled <- fit$intercepts_scaled
#>     beta <- fit$betas
#>     beta_scaled <- fit$betas_scaled
#>     nonzeros <- lapply(beta, function(b) abs(b) > 0)
#>     coefficients <- beta
#>     names(coefficients) <- paste0("p", seq_along(beta))
#>     passes <- fit$passes
#>     reached_max_passes <- passes >= max_passes
#>     if (any(reached_max_passes)) {
#>         reached_max_passes_where <- which(reached_max_passes)
#>         warning("maximum number of passes reached at steps ", 
#>             paste(reached_max_passes_where, collapse = ", "), 
#>             "!")
#>     }
#>     diagnostics <- if (diagnostics) 
#>         setup_diagnostics(fit)
#>     else NULL
#>     patterns <- if (control$patterns) {
#>         fit$patterns
#>     }
#>     else {
#>         NULL
#>     }
#>     slope_class <- switch(control$family, gaussian = "GaussianSLOPE", 
#>         binomial = "BinomialSLOPE", poisson = "PoissonSLOPE", 
#>         multinomial = "MultinomialSLOPE")
#>     structure(list(intercepts = intercepts, coefficients = coefficients, 
#>         intercepts_scaled = intercepts_scaled, coefficients_scaled = beta_scaled, 
#>         nonzeros = nonzeros, lambda = lambda, alpha = alpha[seq_along(beta)], 
#>         variable_names = control$variable_names, class_names = control$class_names, 
#>         passes = passes, deviance_ratio = drop(fit$deviance_ratio), 
#>         null_deviance = fit$null_deviance, family = control$family, 
#>         diagnostics = diagnostics, patterns = patterns, has_intercept = control$fit_intercept, 
#>         n_observations = NROW(x), n_predictors = NCOL(x), call = ocall), 
#>         class = c(slope_class, "SLOPE"))
#> })(x = structure(c(21, 21, 22.8, 21.4, 18.7, 18.1, 14.3, 24.4, 
#> 22.8, 19.2, 17.8, 16.4, 17.3, 15.2, 10.4, 10.4, 14.7, 32.4, 30.4, 
#> 33.9, 21.5, 15.5, 15.2, 13.3, 19.2, 27.3, 26, 30.4, 15.8, 19.7, 
#> 15, 21.4, 3.9, 3.9, 3.85, 3.08, 3.15, 2.76, 3.21, 3.69, 3.92, 
#> 3.92, 3.92, 3.07, 3.07, 3.07, 2.93, 3, 3.23, 4.08, 4.93, 4.22, 
#> 3.7, 2.76, 3.15, 3.73, 3.08, 4.08, 4.43, 3.77, 4.22, 3.62, 3.54, 
#> 4.11, 2.62, 2.875, 2.32, 3.215, 3.44, 3.46, 3.57, 3.19, 3.15, 
#> 3.44, 3.44, 4.07, 3.73, 3.78, 5.25, 5.424, 5.345, 2.2, 1.615, 
#> 1.835, 2.465, 3.52, 3.435, 3.84, 3.845, 1.935, 2.14, 1.513, 3.17, 
#> 2.77, 3.57, 2.78), dim = c(32L, 3L), dimnames = list(c("Mazda RX4", 
#> "Mazda RX4 Wag", "Datsun 710", "Hornet 4 Drive", "Hornet Sportabout", 
#> "Valiant", "Duster 360", "Merc 240D", "Merc 230", "Merc 280", 
#> "Merc 280C", "Merc 450SE", "Merc 450SL", "Merc 450SLC", "Cadillac Fleetwood", 
#> "Lincoln Continental", "Chrysler Imperial", "Fiat 128", "Honda Civic", 
#> "Toyota Corolla", "Toyota Corona", "Dodge Challenger", "AMC Javelin", 
#> "Camaro Z28", "Pontiac Firebird", "Fiat X1-9", "Porsche 914-2", 
#> "Lotus Europa", "Ford Pantera L", "Ferrari Dino", "Maserati Bora", 
#> "Volvo 142E"), c("mpg", "drat", "wt"))), y = structure(c(110, 
#> 110, 93, 110, 175, 105, 245, 62, 95, 123, 123, 180, 180, 180, 
#> 205, 215, 230, 66, 52, 65, 97, 150, 150, 245, 175, 66, 91, 113, 
#> 264, 175, 335, 109), dim = c(32L, 1L)), alpha = 3.48887231288147, 
#>     q = 0.1, gamma = 0)
#> 
#> Path summary:
#>      alpha deviance_ratio n_nonzero
#> 1 3.488872      0.6024373         1

# Or use refit() to refit with different measure
fit <- refit(
  tune,
  subset(mtcars, select = c("mpg", "drat", "wt")),
  mtcars$hp
)
coef(fit)
#> 4 x 1 sparse Matrix of class "dgCMatrix"
#>                
#> [1,] 324.082314
#> [2,]  -8.829731
#> [3,]   .       
#> [4,]   .       
```
