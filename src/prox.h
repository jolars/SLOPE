#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

inline mat prox(const mat& beta, const vec& lambda)
{
  uword p = beta.n_elem;

  // collect sign of beta and work with sorted absolutes
  vec beta_vec = vectorise(beta);
  vec beta_sign = sign(beta_vec);
  beta_vec = abs(beta_vec);
  uvec beta_order = sort_index(beta_vec, "descend");
  beta_vec = (beta_vec(beta_order)).eval();

  vec s(p);
  vec w(p);

  uvec idx_i(p);
  uvec idx_j(p);

  uword k = 0;

  for (uword i = 0; i < p; i++) {
    idx_i(k) = i;
    idx_j(k) = i;
    s(k)     = beta_vec(i) - lambda(i);
    w(k)     = s(k);

    while ((k > 0) && (w(k - 1) <= w(k))) {
      k--;
      idx_j(k)  = i;
      s(k)     += s(k + 1);
      w(k)      = s(k) / (i - idx_i(k) + 1.0);
    }
    k++;
  }

  for (uword j = 0; j < k; j++) {
    double d = std::max(w(j), 0.0);
    for (uword i = idx_i(j); i <= idx_j(j); i++) {
      beta_vec(i) = d;
    }
  }

  // reset order
  beta_vec(beta_order) = beta_vec;

  beta_vec %= beta_sign;

  // reset sign and return
  return reshape(beta_vec, size(beta));
}

