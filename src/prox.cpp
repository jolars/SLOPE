#include "prox.h"

// Stack-based algorithm (Algorithm 4 in Bogdan et al. 2015)
void
prox_stack(arma::vec& x, const arma::vec& lambda)
{
  using namespace arma;

  uword p = x.n_elem;

  vec s(p);
  vec w(p);

  uvec idx_i(p);
  uvec idx_j(p);

  uword k = 0;

  for (uword i = 0; i < p; i++) {
    idx_i(k) = i;
    idx_j(k) = i;
    s(k)     = x(i) - lambda(i);
    w(k)     = s(k);

    while ((k > 0) && (w(k - 1) <= w(k))) {
      k--;
      idx_j(k) = i;
      s(k) += s(k + 1);
      w(k) = s(k) / (i - idx_i(k) + 1.0);
    }
    k++;
  }

  for (uword j = 0; j < k; j++) {
    double d = std::max(w(j), 0.0);
    for (uword i = idx_i(j); i <= idx_j(j); i++) {
      x(i) = d;
    }
  }
}

void
prox_pava(arma::vec& y, const arma::vec& lambda)
{
  using namespace arma;

  uword n = y.n_elem;

  vec yc(n + 1);
  double tmp{ 0 };

  double slope{ 0 };

  yc(0)      = 0;
  yc.tail(n) = cumsum(y - lambda);

  uword known = 0;
  uword ip    = 0;

  do {
    slope = -datum::inf;

    for (uword i = known + 1; i <= n; ++i) {
      tmp = (yc(i) - yc(known)) / (i - known);

      if (tmp > slope) {
        slope = tmp;
        ip    = i;
      }
    }

    for (uword i = known; i < ip; ++i) {
      y(i) = (yc(ip) - yc(known)) / (ip - known);
    }

  } while ((known = ip) < n);

  y = clamp(y, 0.0, datum::inf);
}

arma::mat
prox(const arma::mat& beta,
     const arma::vec& lambda,
     const ProxMethod prox_method)
{
  using namespace arma;

  // collect sign of beta and work with sorted absolutes
  vec beta_vec  = vectorise(beta);
  vec beta_sign = sign(beta_vec);
  beta_vec      = abs(beta_vec);
  uvec ord      = sort_index(beta_vec, "descend");
  beta_vec      = (beta_vec(ord)).eval();

  switch (prox_method) {
    case ProxMethod::stack:
      prox_stack(beta_vec, lambda);
      break;
    case ProxMethod::pava:
      prox_pava(beta_vec, lambda);
      break;
  }

  // reset order
  beta_vec(ord) = beta_vec;

  beta_vec %= beta_sign;

  // reset sign and return
  return reshape(beta_vec, size(beta));
}
