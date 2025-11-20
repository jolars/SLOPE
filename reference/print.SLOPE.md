# Print results from SLOPE fit

Print results from SLOPE fit

## Usage

``` r
# S3 method for class 'SLOPE'
print(x, ...)

# S3 method for class 'TrainedSLOPE'
print(x, ...)
```

## Arguments

- x:

  an object of class `'SLOPE'` or `'TrainedSLOPE'`

- ...:

  other arguments passed to
  [`print()`](https://rdrr.io/r/base/print.html)

## Value

Prints output on the screen

## See also

[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md),
`print.SLOPE()`

Other SLOPE-methods:
[`coef.SLOPE()`](https://jolars.github.io/SLOPE/reference/coef.SLOPE.md),
[`deviance.SLOPE()`](https://jolars.github.io/SLOPE/reference/deviance.SLOPE.md),
[`plot.SLOPE()`](https://jolars.github.io/SLOPE/reference/plot.SLOPE.md),
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md),
[`score()`](https://jolars.github.io/SLOPE/reference/score.md),
[`summary.SLOPE()`](https://jolars.github.io/SLOPE/reference/summary.SLOPE.md)

## Examples

``` r
fit <- SLOPE(wine$x, wine$y, family = "multinomial")
print(fit, digits = 1)
#> 
#> Call:
#> SLOPE(x = wine$x, y = wine$y, family = "multinomial")
#> 
#> Path summary:
#>    alpha deviance_ratio n_nonzero
#> 1  1e-01           0.00         0
#> 2  1e-01           0.08         4
#> 3  1e-01           0.16         4
#> 4  1e-01           0.22         4
#> 5  9e-02           0.28         4
#> 6  8e-02           0.33         4
#> 7  8e-02           0.38         4
#> 8  7e-02           0.42         5
#> 9  6e-02           0.46         5
#> 10 6e-02           0.50        10
#> 11 5e-02           0.54        10
#> 12 5e-02           0.57        13
#> 13 4e-02           0.61        13
#> 14 4e-02           0.64        14
#> 15 4e-02           0.67        14
#> 16 3e-02           0.70        14
#> 17 3e-02           0.72        14
#> 18 3e-02           0.74        14
#> 19 3e-02           0.76        14
#> 20 2e-02           0.78        14
#> 21 2e-02           0.80        14
#> 22 2e-02           0.81        15
#> 23 2e-02           0.83        16
#> 24 2e-02           0.84        16
#> 25 1e-02           0.85        16
#> 26 1e-02           0.86        16
#> 27 1e-02           0.87        16
#> 28 1e-02           0.88        16
#> 29 1e-02           0.89        16
#> 30 9e-03           0.90        16
#> 31 8e-03           0.90        16
#> 32 8e-03           0.91        16
#> 33 7e-03           0.92        16
#> 34 6e-03           0.92        16
#> 35 6e-03           0.93        17
#> 36 5e-03           0.93        17
#> 37 5e-03           0.94        17
#> 38 4e-03           0.94        18
#> 39 4e-03           0.94        18
#> 40 4e-03           0.95        19
#> 41 3e-03           0.95        19
#> 42 3e-03           0.96        19
#> 43 3e-03           0.96        19
#> 44 2e-03           0.96        18
#> 45 2e-03           0.96        18
#> 46 2e-03           0.97        18
#> 47 2e-03           0.97        18
#> 48 2e-03           0.97        18
#> 49 2e-03           0.97        18
#> 50 1e-03           0.98        18
#> 51 1e-03           0.98        18
#> 52 1e-03           0.98        18
#> 53 1e-03           0.98        18
#> 54 1e-03           0.98        18
#> 55 9e-04           0.98        18
#> 56 8e-04           0.98        18
#> 57 7e-04           0.99        18
#> 58 7e-04           0.99        18
#> 59 6e-04           0.99        18
#> 60 6e-04           0.99        18
#> 61 5e-04           0.99        18
#> 62 5e-04           0.99        18
#> 63 4e-04           0.99        18
#> 64 4e-04           0.99        18
#> 65 3e-04           0.99        18
#> 66 3e-04           0.99        19
#> 67 3e-04           0.99        19
#> 68 3e-04           0.99        19
#> 69 2e-04           0.99        19
#> 70 2e-04           1.00        19
#> 71 2e-04           1.00        19
#> 72 2e-04           1.00        19
#> 73 2e-04           1.00        19
#> 74 2e-04           1.00        19
#> 75 1e-04           1.00        19
#> 76 1e-04           1.00        19
#> 77 1e-04           1.00        19
#> 78 1e-04           1.00        19
#> 79 1e-04           1.00        19
#> 80 9e-05           1.00        19
#> 81 8e-05           1.00        19
#> 82 7e-05           1.00        19
#> 83 7e-05           1.00        19
#> 84 6e-05           1.00        19
#> 85 5e-05           1.00        19
#> 86 5e-05           1.00        19
#> 87 5e-05           1.00        19
#> 88 4e-05           1.00        19
```
