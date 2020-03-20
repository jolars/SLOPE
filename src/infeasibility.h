#pragma once

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

inline double infeasibility(const mat& gradient, const vec& lambda)
{
  vec infeas = cumsum(sort(abs(vectorise(gradient)), "descending") - lambda);
  return std::max(infeas.max(), 0.0);
}
