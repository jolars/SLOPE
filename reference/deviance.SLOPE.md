# Model Deviance

Model Deviance

## Usage

``` r
# S3 method for class 'SLOPE'
deviance(object, ...)
```

## Arguments

- object:

  an object of class `'SLOPE'`.

- ...:

  ignored

## Value

For Gaussian models this is twice the residual sums of squares. For all
other models, two times the negative loglikelihood is returned.

## See also

[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md)

Other SLOPE-methods:
[`coef.SLOPE()`](https://jolars.github.io/SLOPE/reference/coef.SLOPE.md),
[`plot.SLOPE()`](https://jolars.github.io/SLOPE/reference/plot.SLOPE.md),
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md),
[`print.SLOPE()`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md),
[`score()`](https://jolars.github.io/SLOPE/reference/score.md),
[`summary.SLOPE()`](https://jolars.github.io/SLOPE/reference/summary.SLOPE.md)

## Examples

``` r
fit <- SLOPE(heart$x, heart$y, family = "binomial")
deviance(fit)
#>  [1] 1.3739232 1.2937213 1.2250313 1.1670221 1.1161942 1.0717583 1.0275800
#>  [8] 0.9889880 0.9549637 0.9247657 0.8987407 0.8756350 0.8552501 0.8352811
#> [15] 0.8140560 0.7943376 0.7767137 0.7611074 0.7457175 0.7323982 0.7197023
#> [22] 0.7086136 0.6984955 0.6897198 0.6822010 0.6754331 0.6691962 0.6637329
#> [29] 0.6587840 0.6545040 0.6506893 0.6473829 0.6445254 0.6420727 0.6400310
#> [36] 0.6381618 0.6365580 0.6350098 0.6336419 0.6324513 0.6314431 0.6305687
#> [43] 0.6298194 0.6291763 0.6286303 0.6281727 0.6277802 0.6274314 0.6271389
#> [50] 0.6268910 0.6266718 0.6264895 0.6263366 0.6262057 0.6260965 0.6260051
#> [57] 0.6259285 0.6258649 0.6258109 0.6257664 0.6257289 0.6256976 0.6256716
#> [64] 0.6256498 0.6256316 0.6256165 0.6256039 0.6255934 0.6255846 0.6255774
#> [71] 0.6255713
```
