#pragma once

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

uvec
kktCheck(mat gradient,
         mat beta,
         const vec& lambda,
         const double tol,
         const bool intercept)
{
  if (intercept) {
    gradient.shed_row(0);
    beta.shed_row(0);
  }

  if (beta.n_rows == 0) {
    uvec tmp;
    return tmp;
  }

  uvec nonzeros = find(beta != 0);
  uvec ord = sort_index(abs(gradient), "descend");
  vec abs_gradient_sorted = abs(gradient(ord));

  double rh = std::max(std::sqrt(datum::eps), tol * lambda(0));

  uvec tmp = cumsum(abs_gradient_sorted - lambda) > rh;
  tmp(ord) = tmp;
  tmp(nonzeros).zeros();

  umat out_mat = reshape(tmp, size(gradient));

  uvec out = find(any(out_mat, 1));

  if (intercept)
    out += 1;

  return out;
}
