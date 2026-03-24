# Compute One of Several Loss Metrics on a New Data Set

This function is a unified interface to return various types of loss for
a model fit with
[`SLOPE()`](https://jolars.github.io/SLOPE/dev/reference/SLOPE.md).

## Usage

``` r
score(object, x, y, measure)

# S3 method for class 'GaussianSLOPE'
score(object, x, y, measure = c("mse", "mae"))

# S3 method for class 'BinomialSLOPE'
score(object, x, y, measure = c("mse", "mae", "deviance", "misclass", "auc"))

# S3 method for class 'MultinomialSLOPE'
score(object, x, y, measure = c("mse", "mae", "deviance", "misclass"))

# S3 method for class 'PoissonSLOPE'
score(object, x, y, measure = c("mse", "mae"))
```

## Arguments

- object:

  an object of class `"SLOPE"`

- x:

  feature matrix

- y:

  response

- measure:

  type of target measure. `"mse"` returns mean squared error. `"mae"`
  returns mean absolute error, `"misclass"` returns misclassification
  rate, and `"auc"` returns area under the ROC curve.

## Value

The measure along the regularization path depending on the value in
`measure`.#'

## See also

[`SLOPE()`](https://jolars.github.io/SLOPE/dev/reference/SLOPE.md),
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/dev/reference/predict.SLOPE.md)

Other SLOPE-methods:
[`coef.SLOPE()`](https://jolars.github.io/SLOPE/dev/reference/coef.SLOPE.md),
[`deviance.SLOPE()`](https://jolars.github.io/SLOPE/dev/reference/deviance.SLOPE.md),
[`plot.SLOPE()`](https://jolars.github.io/SLOPE/dev/reference/plot.SLOPE.md),
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/dev/reference/predict.SLOPE.md),
[`print.SLOPE()`](https://jolars.github.io/SLOPE/dev/reference/print.SLOPE.md),
[`summary.SLOPE()`](https://jolars.github.io/SLOPE/dev/reference/summary.SLOPE.md)

## Examples

``` r
x <- subset(infert, select = c("induced", "age", "pooled.stratum"))
y <- infert$case

fit <- SLOPE(x, y, family = "binomial")
score(fit, x, y, measure = "auc")
#>  [1] 0.4546915 0.4673238 0.5060241 0.5118656 0.4935378 0.4780577 0.4980650
#>  [8] 0.5120117 0.4925885 0.5171961
```
