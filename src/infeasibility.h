#pragma once

#include <RcppArmadillo.h>

inline double
infeasibility(const arma::mat gradient, const arma::vec& lambda)
{
  using namespace arma;

  vec infeas = cumsum(sort(abs(vectorise(gradient)), "descending") - lambda);

  return std::max(infeas.max(), 0.0);
}
