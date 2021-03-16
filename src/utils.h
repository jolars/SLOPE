#pragma once

#include <RcppArmadillo.h>

using namespace arma;

mat matrixSubset(const mat& x, const uvec& active_set)
{
  return x.cols(active_set);
}

sp_mat matrixSubset(const sp_mat& x, const uvec& active_set)
{
  const uword p = active_set.n_elem;
  const uword n = x.n_rows;

  sp_mat x_subset(n, p);

  for (uword j = 0; j < p; ++j) {
    uword k = active_set(j);
    x_subset.col(j) = x.col(k);
  }

  return x_subset;
}

inline uvec setUnion(const uvec& a, const uvec& b)
{
  std::vector<unsigned> out;
  std::set_union(a.begin(), a.end(),
                 b.begin(), b.end(),
                 std::back_inserter(out));

  return conv_to<uvec>::from(out);
}


inline uvec setDiff(uvec& a, uvec& b)
{
  std::vector<unsigned> out;
  std::set_difference(a.begin(), a.end(),
                      b.begin(), b.end(),
                      std::back_inserter(out));

  return conv_to<uvec>::from(out);
}

inline bool isSparse(SEXP x)
{
  bool is_sparse = false;

  if (Rf_isS4(x))
    if (Rf_inherits(x, "dgCMatrix"))
      is_sparse = true;

  return is_sparse;
}
