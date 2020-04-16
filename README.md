
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SLOPE <a href="https://jolars.github.io/SLOPE/"><img src='man/figures/logo.svg' align="right" height="139" /></a>

<!-- badges: start -->

[![R CMD Check via
{tic}](https://github.com/jolars/SLOPE/workflows/R%20CMD%20Check%20via%20%7Btic%7D/badge.svg?branch=master)](https://github.com/jolars/SLOPE/actions)
[![Coverage
status](https://codecov.io/gh/jolars/SLOPE/branch/master/graph/badge.svg)](https://codecov.io/github/jolars/SLOPE?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/SLOPE)](https://CRAN.R-project.org/package=SLOPE)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

Efficient implementations for Sorted L-One Penalized Estimation (SLOPE):
generalized linear models regularized with the sorted L1-norm. There is
support for ordinary least-squares regression, binomial regression,
multinomial regression, and poisson regression, as well as both dense
and sparse predictor matrices. In addition, the package features
predictor screening rules that enable efficient solutions to
high-dimensional problems.

## Installation

You can install the current stable release from
[CRAN](https://cran.r-project.org/) with

``` r
install.packages("SLOPE")
```

or the development version from [GitHub](https://github.com/) with

``` r
# install.packages("remotes")
remotes::install_github("jolars/SLOPE")
```

## Versioning

SLOPE uses [semantic versioning](http://semver.org).

## Code of conduct

Please note that the ‘SLOPE’ project is released with a [Contributor
Code of Conduct](https://jolars.github.io/SLOPE/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
