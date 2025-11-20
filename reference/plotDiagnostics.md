# Plot Results from Diagnostics Collected During Model Fitting

This function plots various diagnostics collected during the model
fitting resulting from a call to
[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md) *provided
that `diagnostics = TRUE`*.

## Usage

``` r
plotDiagnostics(
  object,
  ind = max(object$diagnostics$penalty),
  xvar = c("time", "iteration")
)
```

## Arguments

- object:

  an object of class `"SLOPE"`.

- ind:

  either "last"

- xvar:

  what to place on the x axis. `iteration` plots each iteration, `time`
  plots the wall-clock time.

## Value

Invisibly returns NULL. The function is called for its side effect of
producing a plot.

## See also

[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md)

## Examples

``` r
x <- SLOPE(abalone$x, abalone$y, diagnostics = TRUE)
plotDiagnostics(x)
```
