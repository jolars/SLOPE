# Obtain coefficients

This function returns coefficients from a model fit by
[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md).

## Usage

``` r
# S3 method for class 'SLOPE'
coef(
  object,
  alpha = NULL,
  exact = FALSE,
  simplify = TRUE,
  intercept = TRUE,
  scale = c("original", "normalized"),
  sigma,
  ...
)
```

## Arguments

- object:

  an object of class `'SLOPE'`.

- alpha:

  penalty parameter for SLOPE models; if `NULL`, the values used in the
  original fit will be used

- exact:

  if `TRUE` and the given parameter values differ from those in the
  original fit, the model will be refit by calling
  [`stats::update()`](https://rdrr.io/r/stats/update.html) on the object
  with the new parameters. If `FALSE`, the predicted values will be
  based on interpolated coefficients from the original penalty path.

- simplify:

  if `TRUE`, [`base::drop()`](https://rdrr.io/r/base/drop.html) will be
  called before returning the coefficients to drop extraneous dimensions

- intercept:

  whether to include the intercept in the output; only applicable when
  `simplify = TRUE` and an intercept has been fit.

- scale:

  whether to return the coefficients in the original scale or in the
  normalized scale.

- sigma:

  deprecated. Please use `alpha` instead.

- ...:

  arguments that are passed on to
  [`stats::update()`](https://rdrr.io/r/stats/update.html) (and
  therefore also to
  [`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md)) if
  `exact = TRUE` and the given penalty is not in `object`

## Value

Coefficients from the model.

## Details

If `exact = FALSE` and `alpha` is not in `object`, then the returned
coefficients will be approximated by linear interpolation. If
coefficients from another type of penalty sequence (with a different
`lambda`) are required, however, please use
[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md) to refit
the model.

## See also

[`predict.SLOPE()`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md),
[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md)

Other SLOPE-methods:
[`deviance.SLOPE()`](https://jolars.github.io/SLOPE/reference/deviance.SLOPE.md),
[`plot.SLOPE()`](https://jolars.github.io/SLOPE/reference/plot.SLOPE.md),
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md),
[`print.SLOPE()`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md),
[`score()`](https://jolars.github.io/SLOPE/reference/score.md)

## Examples

``` r
fit <- SLOPE(mtcars$mpg, mtcars$vs, path_length = 10)
coef(fit)
#> 2 x 8 sparse Matrix of class "dgCMatrix"
#>                                                                        
#> [1,] 0.4375 -0.27721605 -0.53407168 -0.62638081 -0.65955499 -0.67147717
#> [2,] .       0.03557461  0.04835946  0.05295409  0.05460532  0.05519874
#>                            
#> [1,] -0.6757618 -0.67730159
#> [2,]  0.0554120  0.05548865
coef(fit, scale = "normalized")
#> 2 x 8 sparse Matrix of class "dgCMatrix"
#>                                                                        
#> [1,] 0.4375 0.4375000 0.4375000 0.4375000 0.4375000 0.4375000 0.4375000
#> [2,] .      0.2110296 0.2868697 0.3141252 0.3239204 0.3274406 0.3287056
#>               
#> [1,] 0.4375000
#> [2,] 0.3291603
```
