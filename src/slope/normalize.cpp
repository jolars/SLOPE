#include "normalize.h"

namespace slope {

// TODO: Make this work for sparse beta
std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
rescaleCoefficients(const Eigen::VectorXd& beta0,
                    const Eigen::SparseMatrix<double>& beta,
                    const Eigen::VectorXd& x_centers,
                    const Eigen::VectorXd& x_scales,
                    const bool intercept)
{
  int m = beta0.size();

  bool centering = x_centers.size() > 0;
  bool scaling = x_scales.size() > 0;

  Eigen::MatrixXd beta0_out = beta0;
  Eigen::SparseMatrix<double> beta_out = beta;

  if (centering || scaling) {
    for (int k = 0; k < m; ++k) {
      double x_bar_beta_sum = 0.0;

      for (Eigen::SparseMatrix<double>::InnerIterator it(beta_out, k); it;
           ++it) {
        int j = it.row();

        if (scaling) {
          it.valueRef() /= x_scales(j);
          // beta_out(j, k) /= x_scales(j);
        }
        if (centering) {
          x_bar_beta_sum += x_centers(j) * it.valueRef();
        }
      }

      if (intercept) {
        beta0_out(k) -= x_bar_beta_sum;
      }
    }
  }

  return { beta0_out, beta_out };
}

} // namespace slope
