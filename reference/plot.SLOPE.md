# Plot coefficients

Plot the fitted model's regression coefficients along the regularization
path.

## Usage

``` r
# S3 method for class 'SLOPE'
plot(
  x,
  intercept = FALSE,
  x_variable = c("alpha", "deviance_ratio", "step"),
  magnitudes = FALSE,
  add_labels = FALSE,
  ...
)
```

## Arguments

- x:

  an object of class `"SLOPE"`

- intercept:

  whether to plot the intercept

- x_variable:

  what to plot on the x axis. `"alpha"` plots the scaling parameter for
  the sequence, `"deviance_ratio"` plots the fraction of deviance
  explained, and `"step"` plots step number.

- magnitudes:

  whether to plot the magnitudes of the coefficients

- add_labels:

  whether to add labels (numbers) on the right side of the plot for each
  coefficient

- ...:

  arguments passed to
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html)

## Value

Invisibly returns NULL. The function is called for its side effect of
producing a plot.

## See also

[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md),
[`plotDiagnostics()`](https://jolars.github.io/SLOPE/reference/plotDiagnostics.md)

Other SLOPE-methods:
[`coef.SLOPE()`](https://jolars.github.io/SLOPE/reference/coef.SLOPE.md),
[`deviance.SLOPE()`](https://jolars.github.io/SLOPE/reference/deviance.SLOPE.md),
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md),
[`print.SLOPE()`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md),
[`score()`](https://jolars.github.io/SLOPE/reference/score.md)

## Examples

``` r
fit <- SLOPE(heart$x, heart$y)
plot(fit)
```
