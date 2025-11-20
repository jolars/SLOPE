# Compute one of several loss metrics on a new data set

This function is a unified interface to return various types of loss for
a model fit with
[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md).

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

[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md),
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md)

Other SLOPE-methods:
[`coef.SLOPE()`](https://jolars.github.io/SLOPE/reference/coef.SLOPE.md),
[`deviance.SLOPE()`](https://jolars.github.io/SLOPE/reference/deviance.SLOPE.md),
[`plot.SLOPE()`](https://jolars.github.io/SLOPE/reference/plot.SLOPE.md),
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md),
[`print.SLOPE()`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md),
[`summary.SLOPE()`](https://jolars.github.io/SLOPE/reference/summary.SLOPE.md)

## Examples

``` r
x <- subset(infert, select = c("induced", "age", "pooled.stratum"))
y <- infert$case

fit <- SLOPE(x, y, family = "binomial")
score(fit, x, y, measure = "auc")
#>  [1] 0.5074845 0.5135451 0.4833881 0.4775465 0.5272727 0.5028112 0.5171961
#>  [8] 0.5063162 0.5030303 0.5404162
```
