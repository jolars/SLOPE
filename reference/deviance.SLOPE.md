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
#>  [1] 1.3739232 1.2937211 1.2250312 1.1670203 1.1161940 1.0717581 1.0275799
#>  [8] 0.9889855 0.9549606 0.9247655 0.8987372 0.8756304 0.8552467 0.8352811
#> [15] 0.8140547 0.7943373 0.7767137 0.7611073 0.7457167 0.7323981 0.7197023
#> [22] 0.7086137 0.6984954 0.6897196 0.6822010 0.6754330 0.6691962 0.6637329
#> [29] 0.6587840 0.6545040 0.6506893 0.6473829 0.6445253 0.6420726 0.6400309
#> [36] 0.6381618 0.6365580 0.6350088 0.6336404 0.6324504 0.6314427 0.6305679
#> [43] 0.6298183 0.6291754 0.6286295 0.6281719 0.6277796 0.6274311 0.6271384
#> [50] 0.6268903 0.6266716 0.6264891 0.6263361 0.6262053 0.6260962 0.6260048
#> [57] 0.6259283 0.6258649 0.6258110 0.6257663 0.6257288 0.6256976 0.6256715
#> [64] 0.6256497 0.6256316 0.6256164 0.6256038 0.6255934 0.6255846 0.6255773
#> [71] 0.6255713
```
