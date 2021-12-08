#include "screening.h"
#include <RcppArmadillo.h>

arma::uvec
strongSet(const arma::mat& gradient_prev,
          const arma::vec& lambda,
          const arma::vec& lambda_prev,
          const bool intercept)
{
  using namespace arma;

  const uword m = gradient_prev.n_cols;
  const uword p = lambda.n_elem;
  const vec abs_grad = abs(vectorise(gradient_prev.tail_rows(p / m)));
  const uvec ord = sort_index(abs_grad, "descend");
  const vec tmp = abs_grad(ord) + lambda_prev - 2 * lambda;

  uword i = 0;
  uword k = 0;

  double s = 0;

  while (i + k < p) {
    s += tmp(k + i);

    if (s >= 0) {
      k = k + i + 1;
      i = 0;
      s = 0;
    } else {
      i++;
    }
  }

  uvec active_set(p, fill::zeros);
  active_set.head(k).ones();

  // reset order
  active_set(ord) = active_set;

  umat active_set_mat = reshape(active_set, p / m, m);

  uvec out = find(any(active_set_mat, 1));

  if (intercept) {
    // shuffle all orders by one and add an index for the intercept
    out += 1;
    urowvec top = { 0 };
    out.insert_rows(0, top);
  }

  return out;
}
