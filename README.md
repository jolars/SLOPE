
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

SLOPE uses [semantic versioning](https://semver.org).

## Code of conduct

Please note that the ‘SLOPE’ project is released with a [Contributor
Code of Conduct](https://jolars.github.io/SLOPE/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
