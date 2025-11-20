# Summarize SLOPE model

Produces a summary of a fitted SLOPE model, including information about
the regularization path, model family, and fitted values.

## Usage

``` r
# S3 method for class 'SLOPE'
summary(object, ...)
```

## Arguments

- object:

  an object of class `'SLOPE'`, typically from a call to
  [`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md)

- ...:

  other arguments (currently ignored)

## Value

An object of class `'summary_SLOPE'` with the following components:

- call:

  the call that produced the model

- family:

  the model family

- n_obs:

  number of observations

- n_predictors:

  number of predictors

- has_intercept:

  whether an intercept was fit

- path_length:

  number of steps in the regularization path

- alpha_range:

  range of alpha values in the path

- deviance_ratio_range:

  range of deviance ratios in the path

- null_deviance:

  null deviance

- path_summary:

  data frame summarizing the regularization path

## See also

[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md),
[`print.summary_SLOPE()`](https://jolars.github.io/SLOPE/reference/print.summary_SLOPE.md)

Other SLOPE-methods:
[`coef.SLOPE()`](https://jolars.github.io/SLOPE/reference/coef.SLOPE.md),
[`deviance.SLOPE()`](https://jolars.github.io/SLOPE/reference/deviance.SLOPE.md),
[`plot.SLOPE()`](https://jolars.github.io/SLOPE/reference/plot.SLOPE.md),
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md),
[`print.SLOPE()`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md),
[`score()`](https://jolars.github.io/SLOPE/reference/score.md)

## Examples

``` r
fit <- SLOPE(heart$x, heart$y)
summary(fit)
#> 
#> Call:
#> SLOPE(x = heart$x, y = heart$y) 
#> 
#> Family: gaussian 
#> Observations: 270 
#> Predictors: 18 
#> Intercept: Yes 
#> 
#> Regularization path:
#>   Length: 64 steps
#>   Alpha range: 0.000257 to 0.0902 
#>   Deviance ratio range: 0 to 0.563 
#>   Null deviance: 0.247 
#> 
#> Path summary (first and last 5 steps):
#>     alpha deviance_ratio n_nonzero
#>  0.090200         0.0000         0
#>  0.082200         0.0802         5
#>  0.074900         0.1480         5
#>  0.068200         0.2050         5
#>  0.062200         0.2520         6
#>  . . .
#>  0.000373         0.5630        18
#>  0.000340         0.5630        18
#>  0.000309         0.5630        18
#>  0.000282         0.5630        18
#>  0.000257         0.5630        18

# Multinomial example
fit_multi <- SLOPE(wine$x, wine$y, family = "multinomial")
summary(fit_multi)
#> 
#> Call:
#> SLOPE(x = wine$x, y = wine$y, family = "multinomial") 
#> 
#> Family: multinomial 
#> Observations: 178 
#> Predictors: 13 
#> Intercept: Yes 
#> 
#> Regularization path:
#>   Length: 88 steps
#>   Alpha range: 4.11e-05 to 0.135 
#>   Deviance ratio range: 0 to 0.999 
#>   Null deviance: 2.17 
#> 
#> Path summary (first and last 5 steps):
#>      alpha deviance_ratio n_nonzero
#>  0.1350000         0.0000         0
#>  0.1230000         0.0818         4
#>  0.1120000         0.1590         4
#>  0.1020000         0.2240         4
#>  0.0928000         0.2820         4
#>  . . .
#> 0.0000597         0.9990        19
#> 0.0000544         0.9990        19
#> 0.0000495         0.9990        19
#> 0.0000451         0.9990        19
#> 0.0000411         0.9990        19
```
