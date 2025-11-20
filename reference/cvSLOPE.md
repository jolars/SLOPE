# Tune SLOPE with cross-validation

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

- call:

  the call

## See also

[`plot.TrainedSLOPE()`](https://jolars.github.io/SLOPE/reference/plot.TrainedSLOPE.md)

Other model-tuning:
[`plot.TrainedSLOPE()`](https://jolars.github.io/SLOPE/reference/plot.TrainedSLOPE.md),
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
```
