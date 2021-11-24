#pragma once

#include <RcppArmadillo.h>
#include "family.h"
#include "../results.h"

using namespace Rcpp;
using namespace arma;

class Multinomial : public Family
{
public:
  template<typename... Ts>
  Multinomial(Ts... args)
    : Family(std::forward<Ts>(args)...)
  {}

  double primal(const mat& y, const mat& lin_pred)
  {
    // logsumexp bit
    vec lp_max = max(lin_pred, 1);
    vec lse = trunc_log(exp(-lp_max) +
                        sum(trunc_exp(lin_pred.each_col() - lp_max), 1)) +
              lp_max;

    return accu(lse) - accu(y % lin_pred);
  }

  double dual(const mat& y, const mat& lin_pred)
  {
    vec lp_max = max(lin_pred, 1);
    vec lse = trunc_log(exp(-lp_max) +
                        sum(trunc_exp(lin_pred.each_col() - lp_max), 1)) +
              lp_max;

    return accu(lse) - accu(lin_pred % trunc_exp(lin_pred.each_col() - lse));
  }

  mat pseudoGradient(const mat& y, const mat& lin_pred)
  {
    vec lp_max = max(lin_pred, 1);
    vec lse = trunc_log(exp(-lp_max) +
                        sum(trunc_exp(lin_pred.each_col() - lp_max), 1)) +
              lp_max;

    return trunc_exp(lin_pred.each_col() - lse) - y;
  }

  rowvec fitNullModel(const mat& y, const uword n_classes)
  {
    const uword m = y.n_cols;

    rowvec mu = mean(y);
    rowvec log_mu = trunc_log(mu);

    return log_mu - accu(log_mu + trunc_log(1 - accu(mu))) / (m + 1);
  }

  std::string name() { return "multinomial"; }
};
