#include "hybrid_cd.h"

namespace slope {

std::pair<double, double>
computeClusterGradientAndHessian(const Eigen::MatrixXd& x,
                                 const int j,
                                 const std::vector<int>& s,
                                 const Clusters& clusters,
                                 const Eigen::VectorXd& w,
                                 const Eigen::VectorXd& residual,
                                 const Eigen::VectorXd& x_centers,
                                 const Eigen::VectorXd& x_scales,
                                 const JitNormalization jit_normalization)
{
  int n = x.rows();

  Eigen::VectorXd x_s = Eigen::VectorXd::Zero(n);

  auto s_it = s.cbegin();
  auto c_it = clusters.cbegin(j);

  for (; c_it != clusters.cend(j); ++c_it, ++s_it) {
    int k = *c_it;
    double s_k = *s_it;

    switch (jit_normalization) {
      case JitNormalization::Both:
        x_s += x.col(k) * (s_k / x_scales(k));
        x_s.array() -= x_centers(k) * s_k / x_scales(k);
        break;

      case JitNormalization::Center:
        x_s += x.col(k) * s_k;
        x_s.array() -= x_centers(k) * s_k;
        break;

      case JitNormalization::Scale:
        x_s += x.col(k) * (s_k / x_scales(k));
        break;

      case JitNormalization::None:
        x_s += x.col(k) * s_k;
        break;
    }
  }

  double hessian_j = x_s.cwiseAbs2().dot(w) / n;
  double gradient_j = x_s.cwiseProduct(w).dot(residual) / n;

  return { hessian_j, gradient_j };
}

std::pair<double, double>
computeClusterGradientAndHessian(const Eigen::SparseMatrix<double>& x,
                                 const int j,
                                 const std::vector<int>& s,
                                 const Clusters& clusters,
                                 const Eigen::VectorXd& w,
                                 const Eigen::VectorXd& residual,
                                 const Eigen::VectorXd& x_centers,
                                 const Eigen::VectorXd& x_scales,
                                 const JitNormalization jit_normalization)
{
  int n = x.rows();

  Eigen::VectorXd weighted_residual = w.cwiseProduct(residual);

  std::vector<Eigen::Triplet<double>> triplets;

  Eigen::SparseMatrix<double> x_s(n, 1);

  double offset = 0;

  auto s_it = s.cbegin();
  auto c_it = clusters.cbegin(j);

  for (; c_it != clusters.cend(j); ++c_it, ++s_it) {
    int k = *c_it;
    double s_k = *s_it;

    switch (jit_normalization) {
      case JitNormalization::Center:
        offset += x_centers(k) * s_k;
        break;
      case JitNormalization::Both:
        offset += x_centers(k) * s_k / x_scales(k);
        break;
      case JitNormalization::Scale:
        break;
      case JitNormalization::None:
        break;
    }

    double v = 0;

    for (Eigen::SparseMatrix<double>::InnerIterator it(x, k); it; ++it) {
      switch (jit_normalization) {
        case JitNormalization::Center:
          [[fallthrough]];
        case JitNormalization::None:
          v = it.value() * s_k;
          break;

        case JitNormalization::Scale:
          [[fallthrough]];
        case JitNormalization::Both:
          v = it.value() * s_k / x_scales(k);
          break;
      }

      triplets.emplace_back(it.row(), 0, v);
    }
  }

  x_s.setFromTriplets(triplets.begin(), triplets.end());

  double hess = 1;
  double grad = 0;

  assert(x_s.nonZeros() > 0);

  switch (jit_normalization) {
    case JitNormalization::Center:
      [[fallthrough]];
    case JitNormalization::Both:
      hess = (x_s.col(0).cwiseAbs2().dot(w) - 2 * offset * x_s.col(0).dot(w) +
              std::pow(offset, 2) * w.sum()) /
             n;
      grad = x_s.col(0).dot(weighted_residual) / n -
             offset * weighted_residual.sum() / n;
      break;

    case JitNormalization::Scale:
      [[fallthrough]];
    case JitNormalization::None:
      hess = x_s.col(0).cwiseAbs2().dot(w) / n;
      grad = x_s.col(0).dot(weighted_residual) / n;
      break;
  }

  return { hess, grad };
}

} // namespace slope
