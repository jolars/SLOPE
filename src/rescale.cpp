#include <RcppArmadillo.h>

void
rescale(arma::cube& betas,
        const arma::rowvec& x_center,
        const arma::rowvec& x_scale,
        const arma::rowvec& y_center,
        const arma::rowvec& y_scale,
        const bool intercept)
{
  using namespace arma;

  const uword p = betas.n_rows;
  const uword m = betas.n_cols;
  const uword path_length = betas.n_slices;

  cube x_bar_beta_sum(1, m, path_length, fill::zeros);

  for (uword k = 0; k < m; ++k) {
    for (uword j = static_cast<uword>(intercept); j < p; ++j) {
      betas.tube(j, k) *= y_scale(k) / x_scale(j);
      x_bar_beta_sum.tube(0, k) += x_center(j) * betas.tube(j, k);
    }

    if (intercept)
      betas.tube(0, k) =
        betas.tube(0, k) * y_scale(k) + y_center(k) - x_bar_beta_sum.tube(0, k);
  }
}
