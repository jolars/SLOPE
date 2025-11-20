# Compute One of Several Loss Metrics on a New Data Set

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
#>  [1] 0.4885725 0.5062432 0.5164659 0.5351588 0.5065352 0.4810515 0.4934648
#>  [8] 0.5140562 0.4787879 0.5174151
```
