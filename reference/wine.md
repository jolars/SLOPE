# Wine Cultivars

A data set of results from chemical analysis of wines grown in Italy
from three different cultivars.

## Usage

``` r
wine
```

## Format

178 observations from 13 variables represented as a list consisting of a
categorical response vector `y` with three levels: *A*, *B*, and *C*
representing different cultivars of wine as well as `x`: a sparse
feature matrix of class 'dgCMatrix' with the following variables:

- alcohol:

  alcoholic content

- malic:

  malic acid

- ash:

  ash

- alcalinity:

  alcalinity of ash

- magnesium:

  magnemium

- phenols:

  total phenols

- flavanoids:

  flavanoids

- nonflavanoids:

  nonflavanoid phenols

- proanthocyanins:

  proanthocyanins

- color:

  color intensity

- hue:

  hue

- dilution:

  OD280/OD315 of diluted wines

- proline:

  proline

## Source

Dua, D. and Karra Taniskidou, E. (2017). UCI Machine Learning Repository
<http://archive.ics.uci.edu/ml/>. Irvine, CA: University of California,
School of Information and Computer Science.

<https://raw.githubusercontent.com/hadley/rminds/master/1-data/wine.csv>

<https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass.html#wine>

## See also

Other datasets:
[`abalone`](https://jolars.github.io/SLOPE/reference/abalone.md),
[`bodyfat`](https://jolars.github.io/SLOPE/reference/bodyfat.md),
[`heart`](https://jolars.github.io/SLOPE/reference/heart.md),
[`student`](https://jolars.github.io/SLOPE/reference/student.md)
