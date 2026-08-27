/**
 * @file
 * @brief Diagnostics for SLOPE optimization
 */

#pragma once

#include "jit_normalization.h"
#include "losses/loss.h"
#include "math.h"
#include "ols.h"
#include "sorted_l1_norm.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <algorithm>
#include <memory>

namespace slope {

namespace detail {

template<typename T>
Eigen::MatrixXd
quadraticOlsDualPointFromDesign(const T& x,
                                const Eigen::MatrixXd& y,
                                const bool intercept)
{
  Eigen::MatrixXd theta(y.rows(), y.cols());

  for (Eigen::Index k = 0; k < y.cols(); ++k) {
    auto [beta0, beta] = fitOls(x, Eigen::VectorXd(y.col(k)), intercept);
    theta.col(k) = x * beta - y.col(k);
    if (intercept) {
      theta.col(k).array() += beta0;
      theta.col(k).array() -= theta.col(k).mean();
    }
  }

  return theta;
}

template<typename T>
Eigen::MatrixXd
quadraticOlsDualPoint(const Eigen::MatrixBase<T>& x,
                      const Eigen::MatrixXd& y,
                      const Eigen::VectorXd& x_centers,
                      const Eigen::VectorXd& x_scales,
                      const JitNormalization jit_normalization,
                      const bool intercept)
{
  Eigen::MatrixXd x_effective = x;
  const bool center = jit_normalization == JitNormalization::Center ||
                      jit_normalization == JitNormalization::Both;
  const bool scale = jit_normalization == JitNormalization::Scale ||
                     jit_normalization == JitNormalization::Both;

  if (center && !intercept) {
    x_effective.rowwise() -= x_centers.transpose();
  }
  if (scale) {
    x_effective.array().rowwise() /= x_scales.transpose().array();
  }

  return quadraticOlsDualPointFromDesign(x_effective, y, intercept);
}

template<typename T>
Eigen::MatrixXd
quadraticOlsDualPoint(const Eigen::SparseMatrixBase<T>& x,
                      const Eigen::MatrixXd& y,
                      const Eigen::VectorXd& x_centers,
                      const Eigen::VectorXd& x_scales,
                      const JitNormalization jit_normalization,
                      const bool intercept)
{
  const bool center = jit_normalization == JitNormalization::Center ||
                      jit_normalization == JitNormalization::Both;
  const bool scale = jit_normalization == JitNormalization::Scale ||
                     jit_normalization == JitNormalization::Both;

  Eigen::SparseMatrix<double> x_effective = x;
  if (scale) {
    for (int k = 0; k < x_effective.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(x_effective, k); it;
           ++it) {
        it.valueRef() /= x_scales(it.col());
      }
    }
  }

  if (!center || intercept) {
    return quadraticOlsDualPointFromDesign(x_effective, y, intercept);
  }

  Eigen::VectorXd effective_centers = x_centers;
  if (scale) {
    effective_centers.array() /= x_scales.array();
  }

  Eigen::MatrixXd theta(y.rows(), y.cols());
  for (Eigen::Index k = 0; k < y.cols(); ++k) {
    auto fit = fitOls(x_effective, Eigen::VectorXd(y.col(k)), true);
    const Eigen::VectorXd& beta = fit.second;
    theta.col(k) = x_effective * beta - y.col(k);
    theta.col(k).array() -= effective_centers.dot(beta);
  }

  return theta;
}

} // namespace detail

/**
 * @brief Scales a candidate into the SLOPE dual constraint and evaluates it.
 *
 * @tparam MatrixType The type of the design matrix.
 * @param beta Current coefficient vector, used to determine the dual-gradient
 * dimensions.
 * @param theta Candidate dual point satisfying the loss conjugate domain.
 * @param loss Pointer to the loss function object.
 * @param sl1_norm Sorted L1 norm object.
 * @param lambda Vector of penalty parameters.
 * @param x Design matrix.
 * @param y Response matrix.
 * @param x_centers Vector of feature means for centering.
 * @param x_scales Vector of feature scales for normalization.
 * @param jit_normalization Just-in-time normalization settings.
 * @return The dual objective at the scaled feasible point.
 */
template<typename MatrixType>
double
computeDualFromPoint(const Eigen::VectorXd& beta,
                     Eigen::MatrixXd theta,
                     const std::unique_ptr<Loss>& loss,
                     const SortedL1Norm& sl1_norm,
                     const Eigen::ArrayXd& lambda,
                     const MatrixType& x,
                     const Eigen::MatrixXd& y,
                     const Eigen::VectorXd& x_centers,
                     const Eigen::VectorXd& x_scales,
                     const JitNormalization& jit_normalization)
{
  const int n = x.rows();
  Eigen::VectorXd gradient(beta.size());

  updateGradient(gradient,
                 x,
                 theta,
                 x_centers,
                 x_scales,
                 Eigen::VectorXd::Ones(n),
                 jit_normalization);

  const double dual_norm = sl1_norm.dualNorm(gradient, lambda);
  theta.array() /= std::max(1.0, dual_norm);

  return loss->dual(theta, y, Eigen::VectorXd::Ones(n));
}

/**
 * @brief Computes the dual objective function value for SLOPE optimization
 *
 * @tparam MatrixType The type of the design matrix
 * @param beta Current coefficient vector
 * @param eta Current linear predictor
 * @param loss Pointer to the loss function object
 * @param sl1_norm Sorted L1 norm object
 * @param lambda Vector of penalty parameters
 * @param x Design matrix
 * @param y Response matrix
 * @param x_centers Vector of feature means for centering
 * @param x_scales Vector of feature scales for normalization
 * @param jit_normalization Just-in-time normalization settings
 * @param intercept Boolean indicating if intercept is included in the model
 *
 * @return double The computed dual objective value
 *
 * @details This function computes the dual objective value for the SLOPE
 * optimization problem. It handles both cases with and without intercept terms,
 * applying appropriate normalization and gradient computations.
 */
template<typename MatrixType>
double
computeDual(const Eigen::VectorXd& beta,
            const Eigen::MatrixXd& eta,
            const std::unique_ptr<Loss>& loss,
            const SortedL1Norm& sl1_norm,
            const Eigen::ArrayXd& lambda,
            const MatrixType& x,
            const Eigen::MatrixXd& y,
            const Eigen::VectorXd& x_centers,
            const Eigen::VectorXd& x_scales,
            const JitNormalization& jit_normalization,
            const bool intercept)
{
  Eigen::MatrixXd theta = loss->dualPoint(eta, y, intercept);
  return computeDualFromPoint(beta,
                              theta,
                              loss,
                              sl1_norm,
                              lambda,
                              x,
                              y,
                              x_centers,
                              x_scales,
                              jit_normalization);
}

} // namespace slope
