#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

void standardize(mat& x,
                 rowvec& x_center,
                 rowvec& x_scale,
                 bool intercept,
                 bool center,
                 std::string scale)
{
  const uword p = x.n_cols;

  for (uword j = static_cast<uword>(intercept); j < p; ++j) {
    if (center) {
      x_center(j) = mean(x.col(j));
      x.col(j) -= x_center(j);
    }

    if (scale == "l1") {
      x_scale(j) = norm(x.col(j), 1);
    } else if (scale == "l2") {
      x_scale(j) = norm(x.col(j), 2);
    } else if (scale == "sd") {
      x_scale(j) = stddev(x.col(j));
    } else if (scale == "max") {
      x_scale(j) = x.col(j).max();
    }

    // don't scale zero-variance predictors
    x_scale(j) = x_scale(j) == 0.0 ? 1.0 : x_scale(j);

    if (scale != "none") {
      x.col(j) /= x_scale(j);
    }
  }
}

void standardize(sp_mat& x,
                 rowvec& x_center,
                 rowvec& x_scale,
                 bool intercept,
                 bool center,
                 std::string scale)
{
  const uword p = x.n_cols;
  const uword n = x.n_rows;

  for (uword j = static_cast<uword>(intercept); j < p; ++j) {
    if (scale == "l1") {
      x_scale(j) = norm(x.col(j), 1);
    } else if (scale == "l2") {
      x_scale(j) = norm(x.col(j), 2);
    } else if (scale == "sd") {
      double xbar = accu(x.col(j))/n;
      x_scale(j) = norm(x.col(j) - xbar)/std::sqrt(n - 1);
    } else if (scale == "max") {
      x_scale(j) = x.col(j).max();
    }

    // don't scale zero-variance predictors
    x_scale(j) = x_scale(j) == 0.0 ? 1.0 : x_scale(j);

    if (scale != "none") {
      x.col(j) /= x_scale(j);
    }
  }
}
