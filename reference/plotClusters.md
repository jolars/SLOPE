# Plot cluster structure

Note that this function requires the `patterns` argument to be set to
`TRUE` in the call to
[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md). Calling
this function on a `SLOPE` object without patterns will result in an
error.

## Usage

``` r
plotClusters(
  x,
  plot_signs = FALSE,
  color_clusters = TRUE,
  include_zeroes = TRUE,
  show_alpha = FALSE,
  alpha_steps = NULL,
  palette = "viridis",
  ...
)
```

## Arguments

- x:

  an object of class `'SLOPE'`

- plot_signs:

  logical, indicating whether to plot signs of estimated coefficients on
  the plot

- color_clusters:

  logical, indicating whether the clusters should have different colors

- include_zeroes:

  logical, indicating whether zero variables should be plotted. Default
  to TRUE

- show_alpha:

  logical, indicating whether labels with alpha values or steps in the
  path should be plotted.

- alpha_steps:

  a vector of integer alpha steps to plot. If `NULL`, all the steps are
  plotted.

- palette:

  a character string specifying the color palette to use for the
  clusters. This is passed to
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html).

- ...:

  additional arguments passed to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html).

## Value

Invisibly returns NULL. The function is called for its side effect of
producing a plot.

## See also

[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md),
[`graphics::image()`](https://rdrr.io/r/graphics/image.html),
[`graphics::text()`](https://rdrr.io/r/graphics/text.html).

## Examples

``` r
set.seed(10)
X <- matrix(rnorm(10000), ncol = 10)
colnames(X) <- paste0("X", 1:10)
beta <- c(rep(10, 3), rep(-20, 2), rep(20, 2), rep(0, 3))
Y <- X %*% beta + rnorm(1000)
fit <- SLOPE(X, Y, patterns = TRUE)

plotClusters(fit)

plotClusters(fit, alpha_steps = 1:10)
```
