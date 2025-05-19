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
                          const int ind,
                          const Eigen::MatrixXd& w,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          const double s,
                          const JitNormalization jit_normalization,
                          const int n)
{
  double gradient = 0.0;
  double hessian = 0.0;

  int p = x.cols();

  auto [k, j] = std::div(ind, p);

  // TODO: Avoid these copies
  Eigen::VectorXd residual_v = residual.col(k);
  Eigen::VectorXd w_v = w.col(k);

  switch (jit_normalization) {
    case JitNormalization::Both:
      gradient = s *
                 (x.col(j).cwiseProduct(w_v).dot(residual_v) -
                  w_v.dot(residual_v) * x_centers(j)) /
                 (n * x_scales(j));
      hessian =
        (x.col(j).cwiseAbs2().dot(w_v) - 2 * x_centers(j) * x.col(j).dot(w_v) +
         std::pow(x_centers(j), 2) * w_v.sum()) /
        (std::pow(x_scales(j), 2) * n);
      break;

    case JitNormalization::Center:
      gradient = s *
                 (x.col(j).cwiseProduct(w_v).dot(residual_v) -
                  w_v.dot(residual_v) * x_centers(j)) /
                 n;
      hessian =
        (x.col(j).cwiseAbs2().dot(w_v) - 2 * x_centers(j) * x.col(j).dot(w_v) +
         std::pow(x_centers(j), 2) * w_v.sum()) /
        n;
      break;

    case JitNormalization::Scale:
      gradient =
        s * (x.col(j).cwiseProduct(w_v).dot(residual_v)) / (n * x_scales(j));
      hessian = x.col(j).cwiseAbs2().dot(w_v) / (std::pow(x_scales(j), 2) * n);
      break;

    case JitNormalization::None:
      gradient = s * (x.col(j).cwiseProduct(w_v).dot(residual_v)) / n;
      hessian = x.col(j).cwiseAbs2().dot(w_v) / n;
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
 * @param w Weights
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
template<typename T>
std::pair<double, double>
computeClusterGradientAndHessian(const Eigen::MatrixBase<T>& x,
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
template<typename T>
std::pair<double, double>
computeClusterGradientAndHessian(const Eigen::SparseMatrixBase<T>& x,
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

    for (typename T::InnerIterator it(x.derived(), j); it; ++it) {
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
 * @param w Working weights
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
double
coordinateDescent(Eigen::VectorXd& beta0,
                  Eigen::VectorXd& beta,
                  Eigen::MatrixXd& residual,
                  Clusters& clusters,
                  const Eigen::ArrayXd& lambda,
                  const T& x,
                  const Eigen::MatrixXd& w,
                  const Eigen::VectorXd& x_centers,
                  const Eigen::VectorXd& x_scales,
                  const bool intercept,
                  const JitNormalization jit_normalization,
                  const bool update_clusters)
{
  using namespace Eigen;

  const int n = x.rows();
  const int p = x.cols();
  const int m = residual.cols();

  double max_abs_gradient = 0;

  for (int c_ind = 0; c_ind < clusters.n_clusters(); ++c_ind) {
    double c_old = clusters.coeff(c_ind);

    if (c_old == 0) {
      // We do not update the zero cluster because it can be very large, but
      // often does not change.
      continue;
    }

    int cluster_size = clusters.cluster_size(c_ind);
    std::vector<int> s;
    s.reserve(cluster_size);

    for (auto c_it = clusters.cbegin(c_ind); c_it != clusters.cend(c_ind);
         ++c_it) {
      double s_ind = sign(beta(*c_it));
      s.emplace_back(s_ind);
    }

    double hess = 1;
    double grad = 0;
    VectorXd x_s(n);

    if (cluster_size == 1) {
      int ind = *clusters.cbegin(c_ind);
      std::tie(grad, hess) = computeGradientAndHessian(
        x, ind, w, residual, x_centers, x_scales, s[0], jit_normalization, n);
    } else {
      std::tie(hess, grad) =
        computeClusterGradientAndHessian(x,
                                         c_ind,
                                         s,
                                         clusters,
                                         w,
                                         residual,
                                         x_centers,
                                         x_scales,
                                         jit_normalization);
    }

    max_abs_gradient = std::max(max_abs_gradient, std::abs(grad));

    double c_tilde;
    int new_index;

    if (lambda(0) == 0) {
      // No regularization
      c_tilde = c_old - grad / hess;
      new_index = c_ind;
    } else {
      std::tie(c_tilde, new_index) =
        slopeThreshold(c_old - grad / hess, c_ind, lambda / hess, clusters);
    }

    double c_diff = c_old - c_tilde;

    if (c_diff != 0) {
      auto s_it = s.cbegin();
      auto c_it = clusters.cbegin(c_ind);
      for (; c_it != clusters.cend(c_ind); ++c_it, ++s_it) {
        int ind = *c_it;
        auto [k, j] = std::div(ind, p);
        double s_ind = *s_it;

        // Update coefficient
        beta(ind) = c_tilde * s_ind;

        // Update residual
        switch (jit_normalization) {
          case JitNormalization::Both:
            residual.col(k) -= x.col(j) * (s_ind * c_diff / x_scales(j));
            residual.col(k).array() +=
              x_centers(j) * s_ind * c_diff / x_scales(j);
            break;

          case JitNormalization::Center:
            residual.col(k) -= x.col(j) * (s_ind * c_diff);
            residual.col(k).array() += x_centers(j) * s_ind * c_diff;
            break;

          case JitNormalization::Scale:
            residual.col(k) -= x.col(j) * (s_ind * c_diff / x_scales(j));
            break;

          case JitNormalization::None:
            residual.col(k) -= x.col(j) * (s_ind * c_diff);
            break;
        }
      }
    }

    if (update_clusters) {
      clusters.update(c_ind, new_index, std::abs(c_tilde));
    } else {
      clusters.setCoeff(c_ind, std::abs(c_tilde));
    }
  }

  if (intercept) {
    for (int k = 0; k < residual.cols(); ++k) {
      double beta0_update = residual.col(k).dot(w.col(k)) / n;
      residual.col(k).array() -= beta0_update;
      beta0(k) -= beta0_update;
    }
  }

  return max_abs_gradient;
}

} // namespace slope
