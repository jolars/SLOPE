# Interpolate Coefficients

Interpolate Coefficients

## Usage

``` r
interpolateCoefficients(beta, intercepts, interpolation_list)
```

## Arguments

- beta:

  coefficients

- intercepts:

  intercepts

- interpolation_list:

  a list generated from
  [`interpolatePenalty()`](https://jolars.github.io/SLOPE/reference/interpolatePenalty.md)

## Value

A matrix (or list of matrices) with new coefficients based on linearly
interpolating from new and old lambda values.
