/**
 * @file
 * @brief Diagnostics for SLOPE optimization
 */

#pragma once

#include "jit_normalization.h"
#include "losses/loss.h"
#include "sorted_l1_norm.h"
#include <Eigen/Dense>
#include <memory>
#include <numeric>

namespace slope {

/**
 * @brief Computes the dual objective function value for SLOPE optimization
 *
 * @tparam MatrixType The type of the design matrix
 * @param beta Current coefficient vector
 * @param residual Residual matrix (y - prediction)
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
            const Eigen::MatrixXd& residual,
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
  int n = x.rows();
  int pm = beta.size();

  Eigen::VectorXd gradient(pm);

  std::vector<int> full_set(pm);
  std::iota(full_set.begin(), full_set.end(), 0);

  updateGradient(gradient,
                 x,
                 residual,
                 full_set,
                 x_centers,
                 x_scales,
                 Eigen::VectorXd::Ones(n),
                 jit_normalization);

  Eigen::MatrixXd theta = residual;

  // First compute gradient with potential offset for intercept case
  Eigen::VectorXd dual_gradient = gradient;

  // TODO: Can we avoid this copy? Maybe revert offset afterwards or,
  // alternatively, solve intercept until convergence and then no longer
  // need the offset at all.
  if (intercept) {
    Eigen::VectorXd theta_mean = theta.colwise().mean();
    theta.rowwise() -= theta_mean.transpose();

    offsetGradient(dual_gradient,
                   x,
                   theta_mean,
                   full_set,
                   x_centers,
                   x_scales,
                   jit_normalization);
  }

  // Common scaling operation
  double dual_norm = sl1_norm.dualNorm(dual_gradient, lambda);
  theta.array() /= std::max(1.0, dual_norm);

  double dual = loss->dual(theta, y, Eigen::VectorXd::Ones(n));

  return dual;
}

} // namespace slope
