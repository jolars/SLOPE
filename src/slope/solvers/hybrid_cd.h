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
 * Computes the gradient and Hessian for a cluster of variables in coordinate
 * descent.
 *
 * This function handles the case when multiple variables are in the same
 * cluster (have the same coefficient magnitude), calculating the combined
 * gradient and Hessian needed for the coordinate descent update.
 *
 * @param x Input matrix
 * @param j Cluster index
 * @param s Vector of signs for each variable in the cluster
 * @param clusters The cluster information object
 * @param w Vector of weights
 * @param residual Residual vector
 * @param x_centers Vector of feature centers (means)
 * @param x_scales Vector of feature scales (standard deviations)
 * @param jit_normalization Normalization strategy (Both, Center, Scale, or
 * None)
 *
 * @return std::pair<double, double> containing:
 *         - first: Hessian of the loss function for the cluster
 *         - second: gradient of the loss function for the cluster
 */
std::pair<double, double>
computeClusterGradientAndHessian(const Eigen::MatrixXd& x,
                                 const int j,
                                 const std::vector<int>& s,
                                 const Clusters& clusters,
                                 const Eigen::VectorXd& w,
                                 const Eigen::VectorXd& residual,
                                 const Eigen::VectorXd& x_centers,
                                 const Eigen::VectorXd& x_scales,
                                 const JitNormalization jit_normalization);

/**
 * Computes the gradient and Hessian for a cluster of variables in coordinate
 * descent (sparse matrix version).
 *
 * This overloaded version handles sparse input matrices, optimizing the
 * computation for this data structure.
 *
 * @param x Input sparse matrix
 * @param j Cluster index
 * @param s Vector of signs for each variable in the cluster
 * @param clusters The cluster information object
 * @param w Vector of weights
 * @param residual Residual vector
 * @param x_centers Vector of feature centers (means)
 * @param x_scales Vector of feature scales (standard deviations)
 * @param jit_normalization Normalization strategy (Both, Center, Scale, or
 * None)
 *
 * @return std::pair<double, double> containing:
 *         - first: Hessian of the loss function for the cluster
 *         - second: gradient of the loss function for the cluster
 */
std::pair<double, double>
computeClusterGradientAndHessian(const Eigen::SparseMatrix<double>& x,
                                 const int j,
                                 const std::vector<int>& s,
                                 const Clusters& clusters,
                                 const Eigen::VectorXd& w,
                                 const Eigen::VectorXd& residual,
                                 const Eigen::VectorXd& x_centers,
                                 const Eigen::VectorXd& x_scales,
                                 const JitNormalization jit_normalization);

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

    int cluster_size = clusters.cluster_size(j);
    std::vector<int> s;
    s.reserve(cluster_size);

    for (auto c_it = clusters.cbegin(j); c_it != clusters.cend(j); ++c_it) {
      double s_k = sign(beta(*c_it));
      s.emplace_back(s_k);
    }

    double hessian_j = 1;
    double gradient_j = 0;
    VectorXd x_s(n);

    if (cluster_size == 1) {
      int k = *clusters.cbegin(j);
      std::tie(gradient_j, hessian_j) = computeGradientAndHessian(
        x, k, w, residual, x_centers, x_scales, s[0], jit_normalization, n);
    } else {
      std::tie(hessian_j, gradient_j) = computeClusterGradientAndHessian(
        x, j, s, clusters, w, residual, x_centers, x_scales, jit_normalization);
    }

    auto [c_tilde, new_index] = slopeThreshold(
      c_old - gradient_j / hessian_j, j, lambda / hessian_j, clusters);

    double c_diff = c_old - c_tilde;

    if (c_diff != 0) {
      auto s_it = s.cbegin();
      auto c_it = clusters.cbegin(j);
      for (; c_it != clusters.cend(j); ++c_it, ++s_it) {
        int k = *c_it;
        double s_k = *s_it;

        // Update coefficient
        beta(k) = c_tilde * s_k;

        // Update residual
        switch (jit_normalization) {
          case JitNormalization::Both:
            residual -= x.col(k) * (s_k * c_diff / x_scales(k));
            residual.array() += x_centers(k) * s_k * c_diff / x_scales(k);
            break;

          case JitNormalization::Center:
            residual -= x.col(k) * (s_k * c_diff);
            residual.array() += x_centers(k) * s_k * c_diff;
            break;

          case JitNormalization::Scale:
            residual -= x.col(k) * (s_k * c_diff / x_scales(k));
            break;

          case JitNormalization::None:
            residual -= x.col(k) * (s_k * c_diff);
            break;
        }
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

} // namespace slope
