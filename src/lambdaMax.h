#pragma once

#include <RcppArmadillo.h>

template<typename T>
arma::vec
lambdaMax(const T& x,
          const arma::mat& y,
          const arma::rowvec& y_scale,
          const arma::uword n_targets,
          const std::string& family,
          const bool intercept)
{
  using namespace arma;

  const uword p = x.n_cols;
  mat lambda_max(p, n_targets);

  if (family == "binomial") {
    vec y_new = (y + 1.0) / 2.0;

    // standardize
    double y_center = mean(y_new);
    y_new -= y_center;

    lambda_max = x.t() * y_new;

  } else if (family == "multinomial") {

    rowvec y_bar = mean(y);
    rowvec y_std = stddev(y, 1);

    mat y_map = y;

    for (uword k = 0; k < n_targets; ++k) {
      y_map.col(k) -= y_bar(k);
      y_map.col(k) /= y_std(k);
    }

    lambda_max = x.t() * y_map;

    for (uword k = 0; k < n_targets; ++k) {
      lambda_max.col(k) *= y_std(k);
    }

  } else if (family == "poisson") {

    lambda_max = x.t() * (1.0 - y);

  } else {

    lambda_max = x.t() * y;
  }

  if (intercept)
    lambda_max.shed_row(0);

  return abs(vectorise(lambda_max));
}
