# Sorted L1 Proximal Operator

The proximal operator for the Sorted L1 Norm, which is the penalty
function in SLOPE. It solves the problem \$\$ \arg\\\min_x \Big(J(x,
\lambda) + \frac{1}{2} \|\|x - v\|\|\_2^2\Big) \$\$ where \\J(x,
\lambda)\\ is the Sorted L1 Norm.

## Usage

``` r
sortedL1Prox(x, lambda, method)
```

## Source

M. Bogdan, E. van den Berg, Chiara Sabatti, Weijie Su, and Emmanuel J.
Candès, “SLOPE – adaptive variable selection via convex optimization,”
Ann Appl Stat, vol. 9, no. 3, pp. 1103–1140, 2015.

## Arguments

- x:

  A vector. In SLOPE, this is the vector of coefficients.

- lambda:

  A non-negative and decreasing sequence of weights for the Sorted L1
  Norm. Needs to be the same length as `x`.

- method:

  DEPRECATED

## Value

An evaluation of the proximal operator at `x` and `lambda`.
