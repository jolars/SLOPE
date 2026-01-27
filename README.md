
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SLOPE <a href="https://jolars.github.io/SLOPE/"><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![R build
status](https://github.com/jolars/SLOPE/workflows/R-CMD-check/badge.svg)](https://github.com/jolars/SLOPE/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/SLOPE)](https://CRAN.R-project.org/package=SLOPE)
[![Code
coverage](https://codecov.io/gh/jolars/SLOPE/graph/badge.svg)](https://app.codecov.io/gh/jolars/SLOPE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17475845.svg)](https://doi.org/10.5281/zenodo.17475845)
<!-- badges: end -->

Efficient implementations for Sorted L-One Penalized Estimation (SLOPE):
generalized linear models regularized with the sorted L1-norm.

## Features

- Gaussian (quadratic), binomial (logistic), multinomial logistic, and
  Poisson regression
- Sparse and dense input matrices
- Efficient hybrid coordinate descent algorithm
- Predictor (feature) screening rules that speed up fitting in
  high-dimensional settings
- Cross-validation
- Parallelized routines
- Duality-based stopping criteria for robust control of suboptimality

## Installation

You can install the current stable release from
[CRAN](https://cran.r-project.org/) with the following command:

``` r
install.packages("SLOPE")
```

Alternatively, you can install the development version from
[GitHub](https://github.com/) with the following command:

``` r
# install.packages("pak")
pak::pak("jolars/SLOPE")
```

## Getting Started

By default, SLOPE fits a full regularization path to the given data.
Here is an example of fitting a logistic SLOPE model to the built-in
`heart` dataset.

``` r
library(SLOPE)

fit <- SLOPE(heart$x, heart$y, family = "binomial")
```

We can plot the resulting regularization path:

``` r
plot(fit)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

We can also perform cross-validation to select optimal scaling of the
regularization sequence:

``` r
set.seed(18)

cvfit <- cvSLOPE(heart$x, heart$y, family = "binomial")
plot(cvfit)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

## Ecosystem

SLOPE is also available as a

- Python package: [sortedl1](https://pypi.org/project/sortedl1/)
- Julia package: [SLOPE.jl](https://github.com/jolars/SLOPE.jl)
- C++ library: [libslope](https://github.com/jolars/libslope)

## Versioning

SLOPE uses [semantic versioning](https://semver.org).

## Code of conduct

Please note that the ‘SLOPE’ project is released with a [Contributor
Code of Conduct](https://jolars.github.io/SLOPE/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
