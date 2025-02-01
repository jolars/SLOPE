#pragma once

#include "family.h"
#include <RcppArmadillo.h>

class Poisson : public Family
{
public:
  template<typename... Ts>
  Poisson(Ts... args)
    : Family(std::forward<Ts>(args)...)
  {
  }

  double primal(const arma::mat& y, const arma::mat& lin_pred)
  {
    using namespace arma;

    return -accu(y % lin_pred - trunc_exp(lin_pred) - lgamma(y + 1));
  }

  double dual(const arma::mat& y, const arma::mat& lin_pred)
  {
    using namespace arma;

    return -accu(trunc_exp(lin_pred) % (lin_pred - 1) - lgamma(y + 1));
  }

  arma::mat partialGradient(const arma::mat& y, const arma::mat& lin_pred)
  {
    return arma::trunc_exp(lin_pred) - y;
  }

  arma::rowvec fitNullModel(const arma::mat& y, const arma::uword n_classes)
  {
    return arma::trunc_log(arma::mean(y));
  }

  std::string name() { return "poisson"; }
};
