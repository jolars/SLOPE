#pragma once

#include <RcppArmadillo.h>

inline arma::mat
matrixSubset(const arma::mat& x, const arma::uvec& active_set)
{
  return x.cols(active_set);
}

inline arma::sp_mat
matrixSubset(const arma::sp_mat& x, const arma::uvec& active_set)
{
  using namespace arma;

  const uword p = active_set.n_elem;
  const uword n = x.n_rows;

  sp_mat x_subset(n, p);

  for (uword j = 0; j < p; ++j) {
    uword k = active_set(j);
    x_subset.col(j) = x.col(k);
  }

  return x_subset;
}

inline arma::uvec
setUnion(const arma::uvec& a, const arma::uvec& b)
{
  std::vector<unsigned> out;
  std::set_union(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  return arma::conv_to<arma::uvec>::from(out);
}

inline arma::uvec
setDiff(arma::uvec& a, arma::uvec& b)
{
  std::vector<unsigned> out;
  std::set_difference(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  return arma::conv_to<arma::uvec>::from(out);
}

inline bool
isSparse(SEXP x)
{
  bool is_sparse = false;

  if (Rf_isS4(x))
    if (Rf_inherits(x, "dgCMatrix"))
      is_sparse = true;

  return is_sparse;
}
