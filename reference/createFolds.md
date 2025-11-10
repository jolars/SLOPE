# Create cross-validation folds

Internal function that creates fold assignments for k-fold
cross-validation, with support for repeated cross-validation.

## Usage

``` r
createFolds(n, n_folds, n_repeats = 1)
```

## Arguments

- n:

  Integer. Number of observations to split into folds.

- n_folds:

  Integer. Number of folds to create.

- n_repeats:

  Integer. Number of times to repeat the cross-validation process with
  different fold assignments. Default is 1.

## Value

A list of length `n_repeats`. Each element contains a list of `n_folds`
integer vectors representing the indices (0-based) of observations in
each fold.
