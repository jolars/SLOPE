#pragma once

#include <RcppArmadillo.h>
#include "family.h"
#include "../results.h"

using namespace Rcpp;
using namespace arma;

class Poisson : public Family
{
public:
  template<typename... Ts>
  Poisson(Ts... args)
    : Family(std::forward<Ts>(args)...)
  {}

  double primal(const mat& y, const mat& lin_pred)
  {
    return -accu(y % lin_pred - trunc_exp(lin_pred) - lgamma(y + 1));
  }

  double dual(const mat& y, const mat& lin_pred)
  {
    return -accu(trunc_exp(lin_pred) % (lin_pred - 1) - lgamma(y + 1));
  }

  mat pseudoGradient(const mat& y, const mat& lin_pred)
  {
    return trunc_exp(lin_pred) - y;
  }

  rowvec fitNullModel(const mat& y, const uword n_classes)
  {
    return trunc_log(mean(y));
  }

  std::string name() { return "poisson"; }
};
