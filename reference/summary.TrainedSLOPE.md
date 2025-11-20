# Summarize TrainedSLOPE Model

Produces a summary of a trained SLOPE model from cross-validation,
including information about the optimal parameters and performance
metrics.

## Usage

``` r
# S3 method for class 'TrainedSLOPE'
summary(object, ...)
```

## Arguments

- object:

  an object of class `'TrainedSLOPE'`, typically from a call to
  [`cvSLOPE()`](https://jolars.github.io/SLOPE/reference/cvSLOPE.md) or
  [`trainSLOPE()`](https://jolars.github.io/SLOPE/reference/trainSLOPE.md)

- ...:

  other arguments (currently ignored)

## Value

An object of class `'summary_TrainedSLOPE'` with the following
components:

- call:

  the call that produced the model

- measure:

  the performance measure(s) used

- optima:

  optimal parameter values and corresponding performance

- n_folds:

  number of cross-validation folds

- n_repeats:

  number of cross-validation repeats

- n_models:

  total number of models evaluated

## See also

[`cvSLOPE()`](https://jolars.github.io/SLOPE/reference/cvSLOPE.md),
[`trainSLOPE()`](https://jolars.github.io/SLOPE/reference/trainSLOPE.md),
[`print.summary_TrainedSLOPE()`](https://jolars.github.io/SLOPE/reference/print.summary_TrainedSLOPE.md)

Other model-tuning:
[`cvSLOPE()`](https://jolars.github.io/SLOPE/reference/cvSLOPE.md),
[`plot.TrainedSLOPE()`](https://jolars.github.io/SLOPE/reference/plot.TrainedSLOPE.md),
[`trainSLOPE()`](https://jolars.github.io/SLOPE/reference/trainSLOPE.md)

## Examples

``` r
tune <- cvSLOPE(
  subset(mtcars, select = c("mpg", "drat", "wt")),
  mtcars$hp,
  q = c(0.1, 0.2),
  n_folds = 5
)
summary(tune)
#> 
#> Call:
#> cvSLOPE(x = subset(mtcars, select = c("mpg", "drat", "wt")), 
#>     y = mtcars$hp, q = c(0.1, 0.2), n_folds = 5) 
#> 
#> Cross-validation:
#>   Folds: 5 
#>   Repeats: 1 
#>   Models evaluated: 133 
#> 
#> Performance measure: Mean Squared Error 
#> 
#> Optimal parameters:
#>    q gamma alpha measure mean  se   lo   hi
#>  0.1     0   2.9     mse 2080 788 -103 4270
```
