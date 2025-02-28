#pragma once

#include "normalize.h"
#include "slope/losses/loss.h"
#include "slope/sorted_l1_norm.h"
#include <Eigen/Dense>
#include <memory>
#include <numeric>

namespace slope {

template<typename MatrixType>
double
computeDual(const Eigen::MatrixXd& beta,
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
  int p = beta.rows();
  int m = beta.cols();

  Eigen::MatrixXd gradient(p, m);

  std::vector<int> full_set(p);
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
  Eigen::MatrixXd dual_gradient = gradient;

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
  double dual_norm = sl1_norm.dualNorm(dual_gradient.reshaped(), lambda);
  theta.array() /= std::max(1.0, dual_norm);

  double dual = loss->dual(theta, y, Eigen::VectorXd::Ones(n));

  return dual;
}

} // namespace slope
