/**
 * @file
 * @brief An implementation of the coordinate descent step in the hybrid
 * algorithm for solving SLOPE.
 */

#pragma once

#include "../clusters.h"
#include "../eigen_compat.h"
#include "../math.h"
#include "slope_threshold.h"
#include <Eigen/Core>
#include <cassert>
#include <random>
#include <vector>

namespace slope {

namespace detail {

class SparseClusterWorkspace
{
public:
  SparseClusterWorkspace() = default;

  SparseClusterWorkspace(const int n, const int m) { resize(n, m); }

  void resize(const int n, const int m)
  {
    if (values.rows() != n || values.cols() != m) {
      values = Eigen::MatrixXd::Zero(n, m);
      active.assign(values.size(), false);
      touched.clear();
    }
  }

  void add(const int row, const int col, const double value)
  {
    const Eigen::Index index = row + values.rows() * col;
    if (!active[index]) {
      active[index] = true;
      touched.emplace_back(index);
    }
    values.data()[index] += value;
  }

  void clear()
  {
    for (const Eigen::Index index : touched) {
      values.data()[index] = 0.0;
      active[index] = false;
    }
    touched.clear();
  }

  Eigen::MatrixXd values;
  std::vector<unsigned char> active;
  std::vector<Eigen::Index> touched;
};

} // namespace detail

/**
 * Computes the gradient and Hessian for coordinate descent optimization with
 * different normalization strategies.
 *
 * @tparam T Matrix type (expected to support col() operations like Eigen
 * matrices)
 *
 * @param x Input matrix
 * @param ind Column index to compute derivatives for
 * @param w Vector of weights
 * @param residual Residual vector
 * @param x_centers Vector of feature centers (means)
 * @param x_scales Vector of feature scales (standard deviations)
 * @param s Step size parameter
 * @param jit_normalization Normalization strategy (Both, Center, Scale, or
 * None)
 * @param n Number of samples
 * @param weight_sums Sum of the working weights for each response
 *
 * @return std::pair<double, double> containing:
 *         - first: gradient of the loss function
 *         - second: diagonal Hessian element
 */
template<typename T>
std::pair<double, double>
computeGradientAndHessian(const Eigen::MatrixBase<T>& x,
                          const int ind,
                          const Eigen::MatrixXd& w,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          const double s,
                          const JitNormalization jit_normalization,
                          const int n,
                          const Eigen::VectorXd& weight_sums)
{
  double gradient = 0.0;
  double hessian = 0.0;

  int p = x.cols();

  auto [k, j] = std::div(ind, p);

  // TODO: Benchmark avoiding these copies in the dense path.
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
         std::pow(x_centers(j), 2) * weight_sums(k)) /
        (std::pow(x_scales(j), 2) * n);
      break;

    case JitNormalization::Center:
      gradient = s *
                 (x.col(j).cwiseProduct(w_v).dot(residual_v) -
                  w_v.dot(residual_v) * x_centers(j)) /
                 n;
      hessian =
        (x.col(j).cwiseAbs2().dot(w_v) - 2 * x_centers(j) * x.col(j).dot(w_v) +
         std::pow(x_centers(j), 2) * weight_sums(k)) /
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

template<typename T>
std::pair<double, double>
computeGradientAndHessian(const Eigen::MatrixBase<T>& x,
                          const int ind,
                          const Eigen::MatrixXd& w,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          const double s,
                          const JitNormalization jit_normalization,
                          const int n)
{
  const Eigen::VectorXd weight_sums = w.colwise().sum().transpose();
  return computeGradientAndHessian(x,
                                   ind,
                                   w,
                                   residual,
                                   x_centers,
                                   x_scales,
                                   s,
                                   jit_normalization,
                                   n,
                                   weight_sums);
}

template<typename T>
std::pair<double, double>
computeGradientAndHessian(const Eigen::SparseMatrixBase<T>& x,
                          const int ind,
                          const Eigen::MatrixXd& w,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          const double s,
                          const JitNormalization jit_normalization,
                          const int n,
                          const Eigen::VectorXd& weight_sums)
{
  const int p = x.cols();
  auto [k, j] = std::div(ind, p);

  double weighted_x_residual_sum = 0.0;
  double weighted_x_sum = 0.0;
  double weighted_x_squared_sum = 0.0;

  for (typename T::InnerIterator it(x.derived(), j); it; ++it) {
    const int i = it.row();
    const double value = it.value();
    const double weight = w(i, k);

    weighted_x_residual_sum += value * weight * residual(i, k);
    weighted_x_sum += value * weight;
    weighted_x_squared_sum += value * value * weight;
  }

  const bool center = jit_normalization == JitNormalization::Center ||
                      jit_normalization == JitNormalization::Both;
  const bool scale = jit_normalization == JitNormalization::Scale ||
                     jit_normalization == JitNormalization::Both;
  const double offset = center ? x_centers(j) : 0.0;
  const double feature_scale = scale ? x_scales(j) : 1.0;

  if (center) {
    weighted_x_residual_sum -=
      offset * w.col(k).cwiseProduct(residual.col(k)).sum();
    weighted_x_squared_sum +=
      offset * offset * weight_sums(k) - 2.0 * offset * weighted_x_sum;
  }

  const double gradient = s * weighted_x_residual_sum / (n * feature_scale);
  const double hessian =
    weighted_x_squared_sum / (n * feature_scale * feature_scale);

  return { gradient, hessian };
}

template<typename T>
std::pair<double, double>
computeGradientAndHessian(const Eigen::SparseMatrixBase<T>& x,
                          const int ind,
                          const Eigen::MatrixXd& w,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          const double s,
                          const JitNormalization jit_normalization,
                          const int n)
{
  const Eigen::VectorXd weight_sums = w.colwise().sum().transpose();
  return computeGradientAndHessian(x,
                                   ind,
                                   w,
                                   residual,
                                   x_centers,
                                   x_scales,
                                   s,
                                   jit_normalization,
                                   n,
                                   weight_sums);
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
 * @param c_ind Cluster index
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
 * @param c_ind Cluster index
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
namespace detail {

template<typename T>
std::pair<double, double>
computeSparseClusterGradientAndHessian(const Eigen::SparseMatrixBase<T>& x,
                                       const int c_ind,
                                       const std::vector<int>& s,
                                       const Clusters& clusters,
                                       const Eigen::MatrixXd& w,
                                       const Eigen::MatrixXd& residual,
                                       const Eigen::VectorXd& weight_sums,
                                       const Eigen::VectorXd& x_centers,
                                       const Eigen::VectorXd& x_scales,
                                       const JitNormalization jit_normalization,
                                       SparseClusterWorkspace& workspace)
{
  int n = x.rows();
  int p = x.cols();
  int m = residual.cols();

  workspace.resize(n, m);

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

    double scale = s_ind;
    if (jit_normalization == JitNormalization::Scale ||
        jit_normalization == JitNormalization::Both) {
      scale /= x_scales(j);
    }

    for (typename T::InnerIterator it(x.derived(), j); it; ++it) {
      workspace.add(it.row(), k, it.value() * scale);
    }
  }

  double hess = 0;
  double grad = 0;
  Eigen::ArrayXd weighted_sum = Eigen::ArrayXd::Zero(m);

  for (const Eigen::Index index : workspace.touched) {
    const int k = index / n;
    const int i = index - k * n;
    const double value = workspace.values.data()[index];
    const double weight = w(i, k);

    hess += value * value * weight;
    grad += value * weight * residual(i, k);
    weighted_sum(k) += value * weight;
  }

  if (jit_normalization == JitNormalization::Center ||
      jit_normalization == JitNormalization::Both) {
    for (int k = 0; k < m; ++k) {
      hess += -2 * offset(k) * weighted_sum(k) +
              std::pow(offset(k), 2) * weight_sums(k);
      grad -= offset(k) * w.col(k).cwiseProduct(residual.col(k)).sum();
    }
  }

  workspace.clear();

  return { hess / n, grad / n };
}

template<typename T>
std::pair<double, double>
computeClusterGradientAndHessianWithWorkspace(
  const Eigen::MatrixBase<T>& x,
  const int c_ind,
  const std::vector<int>& s,
  const Clusters& clusters,
  const Eigen::MatrixXd& w,
  const Eigen::MatrixXd& residual,
  const Eigen::VectorXd&,
  const Eigen::VectorXd& x_centers,
  const Eigen::VectorXd& x_scales,
  const JitNormalization jit_normalization,
  SparseClusterWorkspace&)
{
  return computeClusterGradientAndHessian(
    x, c_ind, s, clusters, w, residual, x_centers, x_scales, jit_normalization);
}

template<typename T>
std::pair<double, double>
computeClusterGradientAndHessianWithWorkspace(
  const Eigen::SparseMatrixBase<T>& x,
  const int c_ind,
  const std::vector<int>& s,
  const Clusters& clusters,
  const Eigen::MatrixXd& w,
  const Eigen::MatrixXd& residual,
  const Eigen::VectorXd& weight_sums,
  const Eigen::VectorXd& x_centers,
  const Eigen::VectorXd& x_scales,
  const JitNormalization jit_normalization,
  SparseClusterWorkspace& workspace)
{
  return computeSparseClusterGradientAndHessian(x,
                                                c_ind,
                                                s,
                                                clusters,
                                                w,
                                                residual,
                                                weight_sums,
                                                x_centers,
                                                x_scales,
                                                jit_normalization,
                                                workspace);
}

} // namespace detail

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
  detail::SparseClusterWorkspace workspace(x.rows(), residual.cols());
  const Eigen::VectorXd weight_sums = w.colwise().sum().transpose();
  return detail::computeSparseClusterGradientAndHessian(x,
                                                        c_ind,
                                                        s,
                                                        clusters,
                                                        w,
                                                        residual,
                                                        weight_sums,
                                                        x_centers,
                                                        x_scales,
                                                        jit_normalization,
                                                        workspace);
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
 * @param lambda_cumsum Cumulative sum of the lambda sequence.
 * @param x The design matrix
 * @param w Working weights
 * @param weight_sums Sum of the working weights for each response
 * @param x_centers The center values of the data matrix columns
 * @param x_scales The scale values of the data matrix columns
 * @param intercept Shuold an intervept be fit?
 * @param jit_normalization Type o fJIT normalization.
 * @param rng Random number generator for shuffling indices in permuted CD.
 * @param update_clusters Whether to maintain cluster ordering and membership
 *   after each update. If false, later cluster-coordinate updates may no longer
 *   satisfy their required invariants.
 * @param cd_type Type of coordinate descent to use ("cyclical" or "permuted")
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
                  const Eigen::ArrayXd& lambda_cumsum,
                  const T& x,
                  const Eigen::MatrixXd& w,
                  const Eigen::VectorXd& weight_sums,
                  const Eigen::VectorXd& x_centers,
                  const Eigen::VectorXd& x_scales,
                  const bool intercept,
                  const JitNormalization jit_normalization,
                  const bool update_clusters,
                  std::mt19937& rng,
                  const std::string& cd_type = "cyclical")
{
  using namespace Eigen;

  const int n = x.rows();
  const int p = x.cols();
  const int m = residual.cols();

  double max_abs_gradient = 0;
  detail::SparseClusterWorkspace sparse_workspace;

  // Create a vector of indices to process
  std::vector<int> indices;
  indices.reserve(clusters.size());
  for (int i = 0; i < clusters.size(); ++i) {
    if (clusters.coeff(i) != 0) { // Skip zero cluster
      indices.push_back(i);
    }
  }

  if (cd_type == "permuted") {
    std::shuffle(indices.begin(), indices.end(), rng);
  }

  for (int c_ind : indices) {
    // Skip if index is no longer valid due to cluster updates
    if (c_ind >= clusters.size()) {
      continue;
    }

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
      int ind = *c_it;
      assert(ind >= 0 && ind < beta.size() && "Invalid index in cluster");
      double s_ind = sign(beta(ind));
      s.emplace_back(s_ind);
    }

    double hess = 1;
    double grad = 0;
    VectorXd x_s(n);

    if (cluster_size == 1) {
      int ind = *clusters.cbegin(c_ind);
      std::tie(grad, hess) = computeGradientAndHessian(x,
                                                       ind,
                                                       w,
                                                       residual,
                                                       x_centers,
                                                       x_scales,
                                                       s[0],
                                                       jit_normalization,
                                                       n,
                                                       weight_sums);
    } else {
      std::tie(hess, grad) =
        detail::computeClusterGradientAndHessianWithWorkspace(x,
                                                              c_ind,
                                                              s,
                                                              clusters,
                                                              w,
                                                              residual,
                                                              weight_sums,
                                                              x_centers,
                                                              x_scales,
                                                              jit_normalization,
                                                              sparse_workspace);
    }

    max_abs_gradient = std::max(max_abs_gradient, std::abs(grad));

    double c_tilde;
    int new_index;

    const double gamma = hess * c_old - grad;
    std::tie(c_tilde, new_index) =
      slopeThreshold(gamma, hess, c_ind, lambda_cumsum, clusters);

    assert(c_tilde == 0 || new_index < clusters.size());
    assert(new_index >= 0 && new_index <= clusters.size());

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
