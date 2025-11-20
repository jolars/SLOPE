# Plot results from cross-validation

Plot results from cross-validation

## Usage

``` r
# S3 method for class 'TrainedSLOPE'
plot(
  x,
  plot_min = TRUE,
  ci_alpha = 0.2,
  ci_border = NA,
  ci_col = "salmon",
  plot_args = list(),
  polygon_args = list(),
  lines_args = list(),
  abline_args = list(),
  index = NULL,
  measure,
  ...
)
```

## Arguments

- x:

  an object of class `'TrainedSLOPE'`, typically from a call to
  [`cvSLOPE()`](https://jolars.github.io/SLOPE/reference/cvSLOPE.md)

- plot_min:

  whether to mark the location of the penalty corresponding to the best
  prediction score

- ci_alpha:

  alpha (opacity) for fill in confidence limits

- ci_border:

  color (or flag to turn off and on) the border of the confidence limits

- ci_col:

  color for border of confidence limits

- plot_args:

  list of additional arguments to pass to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), which sets
  up the plot frame

- polygon_args:

  list of additional arguments to pass to
  [`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html),
  which fills the confidence limits

- lines_args:

  list of additional arguments to pass to
  [`graphics::lines()`](https://rdrr.io/r/graphics/lines.html), which
  plots the mean

- abline_args:

  list of additional arguments to pass to
  [`graphics::abline()`](https://rdrr.io/r/graphics/abline.html), which
  plots the minimum

- index:

  an optional index, to plot only one (the index-th) set of the
  parameter combinations.

- measure:

  any of the measures used in the call to
  [`trainSLOPE()`](https://jolars.github.io/SLOPE/reference/trainSLOPE.md).
  If `measure = "auto"` then deviance will be used for binomial and
  multinomial models, whilst mean-squared error will be used for
  Gaussian and Poisson models.

- ...:

  ignored

## Value

A plot for every value of `q` is produced on the current device.

## See also

[`cvSLOPE()`](https://jolars.github.io/SLOPE/reference/cvSLOPE.md)

Other model-tuning:
[`cvSLOPE()`](https://jolars.github.io/SLOPE/reference/cvSLOPE.md),
[`summary.TrainedSLOPE()`](https://jolars.github.io/SLOPE/reference/summary.TrainedSLOPE.md),
[`trainSLOPE()`](https://jolars.github.io/SLOPE/reference/trainSLOPE.md)

## Examples

``` r
# Cross-validation for a SLOPE binomial model
set.seed(123)
tune <- cvSLOPE(
  subset(mtcars, select = c("mpg", "drat", "wt")),
  mtcars$hp,
  q = c(0.1, 0.2),
  n_folds = 10
)
plot(tune, ci_col = "salmon", index = 1)
```
