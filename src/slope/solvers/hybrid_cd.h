/**
 * @file
 * @brief An implementation of the coordinate descent step in the hybrid
 * algorithm for solving SLOPE.
 */

#pragma once

#include "../clusters.h"
#include "../math.h"
#include "slope_threshold.h"
#include <Eigen/Core>
#include <vector>

namespace slope {
namespace solvers {

/**
 * Computes the gradient and Hessian for coordinate descent optimization with
 * different normalization strategies.
 *
 * @tparam T Matrix type (expected to support col() operations like Eigen
 * matrices)
 *
 * @param x Input matrix
 * @param k Column index to compute derivatives for
 * @param w Vector of weights
 * @param residual Residual vector
 * @param x_centers Vector of feature centers (means)
 * @param x_scales Vector of feature scales (standard deviations)
 * @param s Step size parameter
 * @param jit_normalization Normalization strategy (Both, Center, Scale, or
 * None)
 * @param n Number of samples
 *
 * @return std::pair<double, double> containing:
 *         - first: gradient of the loss function
 *         - second: diagonal Hessian element
 *
 * The function handles four different normalization cases:
 * - Both: Applies both centering and scaling
 * - Center: Applies only centering
 * - Scale: Applies only scaling
 * - None: No normalization
 *
 * Each case computes the gradient and Hessian differently based on the
 * normalization strategy.
 */
template<typename T>
std::pair<double, double>
computeGradientAndHessian(const T& x,
                          const int k,
                          const Eigen::VectorXd& w,
                          const Eigen::VectorXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          const double s,
                          const JitNormalization jit_normalization,
                          const int n)
{
  double gradient = 0.0;
  double hessian = 0.0;

  switch (jit_normalization) {
    case JitNormalization::Both:
      gradient = s *
                 (x.col(k).cwiseProduct(w).dot(residual) -
                  w.dot(residual) * x_centers(k)) /
                 (n * x_scales(k));
      hessian =
        (x.col(k).cwiseAbs2().dot(w) - 2 * x_centers(k) * x.col(k).dot(w) +
         std::pow(x_centers(k), 2) * w.sum()) /
        (std::pow(x_scales(k), 2) * n);
      break;

    case JitNormalization::Center:
      gradient = s *
                 (x.col(k).cwiseProduct(w).dot(residual) -
                  w.dot(residual) * x_centers(k)) /
                 n;
      hessian =
        (x.col(k).cwiseAbs2().dot(w) - 2 * x_centers(k) * x.col(k).dot(w) +
         std::pow(x_centers(k), 2) * w.sum()) /
        n;
      break;

    case JitNormalization::Scale:
      gradient =
        s * (x.col(k).cwiseProduct(w).dot(residual)) / (n * x_scales(k));
      hessian = x.col(k).cwiseAbs2().dot(w) / (std::pow(x_scales(k), 2) * n);
      break;

    case JitNormalization::None:
      gradient = s * (x.col(k).cwiseProduct(w).dot(residual)) / n;
      hessian = x.col(k).cwiseAbs2().dot(w) / n;
      break;
  }

  return { gradient, hessian };
}

/**
 * Coordinate Descent Step
 *
 * This function takes a coordinate descent step in the hybrid CD/PGD algorithm
 * for SLOPE.
 *
 * @tparam T The type of the design matrix. This can be either a dense or
 * sparse.
 * @param beta0 The intercept
 * @param beta The coefficients
 * @param residual The residual vector
 * @param clusters The cluster information, stored in a Cluster object.
 * @param lambda Regularization weights
 * @param x The design matrix
 * @param w The weight vector
 * @param x_centers The center values of the data matrix columns
 * @param x_scales The scale values of the data matrix columns
 * @param intercept Shuold an intervept be fit?
 * @param jit_normalization Type o fJIT normalization.
 * @param update_clusters Flag indicating whether to update the cluster
 * information
 *
 * @see Clusters
 * @see SortedL1Norm
 * @see JitNormalization
 */
template<typename T>
void
coordinateDescent(Eigen::VectorXd& beta0,
                  Eigen::VectorXd& beta,
                  Eigen::VectorXd& residual,
                  Clusters& clusters,
                  const Eigen::ArrayXd& lambda,
                  const T& x,
                  const Eigen::VectorXd& w,
                  const Eigen::VectorXd& x_centers,
                  const Eigen::VectorXd& x_scales,
                  const bool intercept,
                  const JitNormalization jit_normalization,
                  const bool update_clusters)
{
  using namespace Eigen;

  const int n = x.rows();

  for (int j = 0; j < clusters.n_clusters(); ++j) {
    double c_old = clusters.coeff(j);

    if (c_old == 0) {
      // We do not update the zero cluster because it can be very large, but
      // often does not change.
      continue;
    }

    std::vector<int> s;
    int cluster_size = clusters.cluster_size(j);
    s.reserve(cluster_size);

    double hessian_j = 1;
    double gradient_j = 0;
    VectorXd x_s(n);

    if (cluster_size == 1) {
      int k = *clusters.cbegin(j);
      double s_k = sign(beta(k));
      s.emplace_back(s_k);

      std::tie(gradient_j, hessian_j) = computeGradientAndHessian(
        x, k, w, residual, x_centers, x_scales, s_k, jit_normalization, n);
    } else {
      // There's no reasonable just-in-time standardization approach for sparse
      // design matrices when there are clusters in the data, so we need to
      // reduce to a dense column vector.
      x_s.setZero();

      for (auto c_it = clusters.cbegin(j); c_it != clusters.cend(j); ++c_it) {
        int k = *c_it;
        double s_k = sign(beta(k));
        s.emplace_back(s_k);

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

      hessian_j = x_s.cwiseAbs2().dot(w) / n;
      gradient_j = x_s.cwiseProduct(w).dot(residual) / n;
    }

    auto [c_tilde, new_index] = slopeThreshold(
      c_old - gradient_j / hessian_j, j, lambda / hessian_j, clusters);

    auto s_it = s.cbegin();
    auto c_it = clusters.cbegin(j);
    for (; c_it != clusters.cend(j); ++c_it, ++s_it) {
      beta(*c_it) = c_tilde * (*s_it);
    }

    double c_diff = c_old - c_tilde;

    if (c_diff != 0) {
      if (cluster_size == 1) {
        int k = *clusters.cbegin(j);

        switch (jit_normalization) {
          case JitNormalization::Both:
            residual -= x.col(k) * (s[0] * c_diff / x_scales(k));
            residual.array() += x_centers(k) * s[0] * c_diff / x_scales(k);
            break;

          case JitNormalization::Center:
            residual -= x.col(k) * (s[0] * c_diff);
            residual.array() += x_centers(k) * s[0] * c_diff;
            break;

          case JitNormalization::Scale:
            residual -= x.col(k) * (s[0] * c_diff / x_scales(k));
            break;

          case JitNormalization::None:
            residual -= x.col(k) * (s[0] * c_diff);
            break;
        }
      } else {
        residual -= x_s * c_diff;
      }
    }

    if (update_clusters) {
      clusters.update(j, new_index, std::abs(c_tilde));
    } else {
      clusters.setCoeff(j, std::abs(c_tilde));
    }
  }

  if (intercept) {
    double beta0_update = residual.dot(w) / n;
    residual.array() -= beta0_update;
    beta0(0) -= beta0_update;
  }
}

} // namespace solvers
} // namespace slope
