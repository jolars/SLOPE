# Train a SLOPE model

This function trains a model fit by
[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md) by tuning
its parameters through cross-validation.

## Usage

``` r
trainSLOPE(
  x,
  y,
  q = 0.2,
  number = 10,
  repeats = 1,
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

  parameter controlling the shape of the lambda sequence, with usage
  varying depending on the type of path used and has no effect is a
  custom `lambda` sequence is used. Must be greater than `1e-6` and
  smaller than 1.

- number:

  number of folds (cross-validation)

- repeats:

  number of repeats for each fold (for repeated *k*-fold cross
  validation)

- measure:

  measure to try to optimize; note that you may supply *multiple* values
  here and that, by default, all the possible measures for the given
  model will be used.

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

  a `data.frame` listing the used metrics and their labels

- model:

  the model fit to the entire data set

- call:

  the call

## Details

Note that by default this method matches all of the available metrics
for the given model family against those provided in the argument
`measure`. Collecting these measures is not particularly demanding
computationally so it is almost always best to leave this argument as it
is and then choose which argument to focus on in the call to
[`plot.TrainedSLOPE()`](https://jolars.github.io/SLOPE/reference/plot.TrainedSLOPE.md).

## Parallel operation

This function uses the **foreach** package to enable parallel operation.
To enable this, simply register a parallel backend using, for instance,
`doParallel::registerDoParallel()` from the **doParallel** package
before running this function.

## See also

Other model-tuning:
[`cvSLOPE()`](https://jolars.github.io/SLOPE/reference/cvSLOPE.md),
[`plot.TrainedSLOPE()`](https://jolars.github.io/SLOPE/reference/plot.TrainedSLOPE.md)

## Examples

``` r
# 8-fold cross-validation repeated 5 times
tune <- trainSLOPE(subset(mtcars, select = c("mpg", "drat", "wt")),
  mtcars$hp,
  q = c(0.1, 0.2),
  number = 8,
  repeats = 5,
  measure = "mse"
)
```
