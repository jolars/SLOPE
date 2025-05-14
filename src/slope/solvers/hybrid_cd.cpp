#include "hybrid_cd.h"

namespace slope {

std::pair<double, double>
computeClusterGradientAndHessian(const Eigen::MatrixXd& x,
                                 const int c_ind,
                                 const std::vector<int>& s,
                                 const Clusters& clusters,
                                 const Eigen::MatrixXd& w,
                                 const Eigen::MatrixXd& residual,
                                 const Eigen::VectorXd& x_centers,
                                 const Eigen::VectorXd& x_scales,
                                 const JitNormalization jit_normalization)
{
  int n = x.rows();
  int p = x.cols();
  int m = residual.cols();

  Eigen::MatrixXd x_s = Eigen::MatrixXd::Zero(n, m);

  auto s_it = s.cbegin();
  auto c_it = clusters.cbegin(c_ind);

  for (; c_it != clusters.cend(c_ind); ++c_it, ++s_it) {
    int ind = *c_it;
    auto [k, j] = std::div(ind, p);
    double s = *s_it;

    switch (jit_normalization) {
      case JitNormalization::Both:
        x_s.col(k) += x.col(j) * (s / x_scales(j));
        x_s.col(k).array() -= x_centers(j) * s / x_scales(j);
        break;

      case JitNormalization::Center:
        x_s.col(k) += x.col(j) * s;
        x_s.col(k).array() -= x_centers(j) * s;
        break;

      case JitNormalization::Scale:
        x_s.col(k) += x.col(j) * (s / x_scales(j));
        break;

      case JitNormalization::None:
        x_s.col(k) += x.col(j) * s;
        break;
    }
  }

  double hess = 0;
  double grad = 0;

  for (int k = 0; k < m; ++k) {
    hess += x_s.col(k).cwiseAbs2().dot(w.col(k)) / n;
    grad += x_s.col(k).cwiseProduct(w.col(k)).dot(residual.col(k)) / n;
  }

  return { hess, grad };
}

std::pair<double, double>
computeClusterGradientAndHessian(const Eigen::SparseMatrix<double>& x,
                                 const int c_ind,
                                 const std::vector<int>& s,
                                 const Clusters& clusters,
                                 const Eigen::MatrixXd& w,
                                 const Eigen::MatrixXd& residual,
                                 const Eigen::VectorXd& x_centers,
                                 const Eigen::VectorXd& x_scales,
                                 const JitNormalization jit_normalization)
{
  int n = x.rows();
  int p = x.cols();
  int m = residual.cols();

  std::vector<Eigen::Triplet<double>> triplets;

  Eigen::SparseMatrix<double> x_s(n, m);
  Eigen::ArrayXd offset = Eigen::ArrayXd::Zero(m);

  auto s_it = s.cbegin();
  auto c_it = clusters.cbegin(c_ind);

  for (; c_it != clusters.cend(c_ind); ++c_it, ++s_it) {
    int ind = *c_it;
    auto [k, j] = std::div(ind, p);
    double s_ind = *s_it;

    switch (jit_normalization) {
      case JitNormalization::Center:
        offset(k) += x_centers(j) * s_ind;
        break;
      case JitNormalization::Both:
        offset(k) += x_centers(j) * s_ind / x_scales(j);
        break;
      case JitNormalization::Scale:
        break;
      case JitNormalization::None:
        break;
    }

    double v = 0;

    for (Eigen::SparseMatrix<double>::InnerIterator it(x, j); it; ++it) {
      switch (jit_normalization) {
        case JitNormalization::Center:
          [[fallthrough]];
        case JitNormalization::None:
          v = it.value() * s_ind;
          break;

        case JitNormalization::Scale:
          [[fallthrough]];
        case JitNormalization::Both:
          v = it.value() * s_ind / x_scales(j);
          break;
      }

      triplets.emplace_back(it.row(), k, v);
    }
  }

  x_s.setFromTriplets(triplets.begin(), triplets.end());
  assert(x_s.nonZeros() > 0);

  double hess = 0;
  double grad = 0;

  for (int k = 0; k < m; ++k) {
    Eigen::VectorXd weighted_residual = w.col(k).cwiseProduct(residual.col(k));

    switch (jit_normalization) {
      case JitNormalization::Center:
        [[fallthrough]];
      case JitNormalization::Both:
        hess += (x_s.col(k).cwiseAbs2().dot(w.col(k)) -
                 2 * offset(k) * x_s.col(k).dot(w.col(k)) +
                 std::pow(offset(k), 2) * w.col(k).sum()) /
                n;
        grad += x_s.col(k).dot(weighted_residual) / n -
                offset(k) * weighted_residual.sum() / n;
        break;

      case JitNormalization::Scale:
        [[fallthrough]];
      case JitNormalization::None:
        hess += x_s.col(k).cwiseAbs2().dot(w.col(k)) / n;
        grad += x_s.col(k).dot(weighted_residual) / n;

        break;
    }
  }

  return { hess, grad };
}

} // namespace slope
