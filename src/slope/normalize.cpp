#include <slope/normalize.h>

namespace slope {

namespace detail {

double
WeightRangeSums::sum(Eigen::Index begin, Eigen::Index end)
{
  if (begin == end) {
    return 0.0;
  }

  const Eigen::Index n = weights.size();
  if (sums.size() == 0) {
    // Dense columns and scales other than SD never need this cache.
    sums.resize(2 * n);
    sums.tail(n) = weights;
    for (Eigen::Index i = n - 1; i > 0; --i) {
      sums(i) = sums(2 * i) + sums(2 * i + 1);
    }
  }

  if (begin == 0 && end == n) {
    return sums(1);
  }

  // Disjoint tree nodes preserve tiny omitted weights without subtracting
  // nearly equal cumulative sums or rescanning the omitted rows.
  double total = 0.0;
  for (begin += n, end += n; begin < end; begin /= 2, end /= 2) {
    if (begin % 2 != 0) {
      total += sums(begin++);
    }
    if (end % 2 != 0) {
      total += sums(--end);
    }
  }
  return total;
}

} // namespace detail

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
