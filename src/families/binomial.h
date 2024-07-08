#pragma once

#include "family.h"
#include <RcppArmadillo.h>

class Binomial : public Family
{
public:
  template<typename... Ts>
  Binomial(Ts... args)
    : Family(std::forward<Ts>(args)...)
  {}

  double primal(const arma::mat& y, const arma::mat& lin_pred)
  {
    using namespace arma;

    return accu(trunc_log(1.0 + trunc_exp(-y % lin_pred)));
  }

  double dual(const arma::mat& y, const arma::mat& lin_pred)
  {
    using namespace arma;

    const vec r = 1.0 / (1.0 + trunc_exp(y % lin_pred));
    return dot(r - 1.0, trunc_log(1.0 - r)) - dot(r, trunc_log(r));
  }

  arma::mat partialGradient(const arma::mat& y, const arma::mat& lin_pred)
  {
    return -y / (1.0 + arma::trunc_exp(y % lin_pred));
  }

  arma::rowvec fitNullModel(const arma::mat& y, const arma::uword n_classes)
  {
    using namespace arma;

    double pmin = 1e-9;
    double pmax = 1 - pmin;

    vec mu = clamp(mean(0.5 * y + 0.5), pmin, pmax);

    return trunc_log(mu / (1.0 - mu));
  }

  std::string name() { return "binomial"; }
};
