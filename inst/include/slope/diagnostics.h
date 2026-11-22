/**
 * @file
 * @brief Diagnostics for SLOPE optimization
 */

#pragma once

#include "jit_normalization.h"
#include "losses/loss.h"
#include "math.h"
#include "observation_weights.h"
#include "ols.h"
#include "sorted_l1_norm.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <algorithm>
#include <memory>

namespace slope {

namespace detail {

template<typename T>
std::pair<double, Eigen::VectorXd>
fitWeightedOls(const Eigen::MatrixBase<T>& x,
               const Eigen::VectorXd& y,
               const bool intercept,
               const Eigen::VectorXd& weights)
{
  const Eigen::VectorXd root_weights = weights.array().sqrt();
  Eigen::MatrixXd design =
    Eigen::MatrixXd::Zero(x.rows(), x.cols() + static_cast<int>(intercept));
  Eigen::VectorXd response = Eigen::VectorXd::Zero(y.size());
  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    if (weights(i) > 0.0) {
      design.row(i).head(x.cols()) = x.row(i) * root_weights(i);
      response(i) = y(i) * root_weights(i);
    }
  }
  if (intercept) {
    design.col(x.cols()) = root_weights;
  }
  auto fit = fitOls(design, response, false);
  if (!intercept) {
    return fit;
  }
  const double beta0 = fit.second(x.cols());
  return { beta0, fit.second.head(x.cols()) };
}

template<typename T>
std::pair<double, Eigen::VectorXd>
fitWeightedOls(const Eigen::SparseMatrixBase<T>& x,
               const Eigen::VectorXd& y,
               const bool intercept,
               const Eigen::VectorXd& weights)
{
  const Eigen::VectorXd root_weights = weights.array().sqrt();
  Eigen::SparseMatrix<double> weighted_x = x;
  for (int k = 0; k < weighted_x.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(weighted_x, k); it;
         ++it) {
      it.valueRef() =
        weights(it.row()) > 0.0 ? it.value() * root_weights(it.row()) : 0.0;
    }
  }
  Eigen::VectorXd response = Eigen::VectorXd::Zero(y.size());
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    if (weights(i) > 0.0) {
      response(i) = y(i) * root_weights(i);
    }
  }
  Eigen::SparseMatrix<double> design(x.rows(),
                                     x.cols() + static_cast<int>(intercept));
  design.leftCols(x.cols()) = weighted_x;
  if (intercept) {
    Eigen::SparseMatrix<double> intercept_column = root_weights.sparseView();
    design.rightCols(1) = intercept_column;
  }
  auto fit = fitOls(design, response, false);
  if (!intercept) {
    return fit;
  }
  const double beta0 = fit.second(x.cols());
  return { beta0, fit.second.head(x.cols()) };
}

template<typename T>
Eigen::MatrixXd
quadraticOlsDualPointFromDesign(
  const T& x,
  const Eigen::MatrixXd& y,
  const bool intercept,
  const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
{
  const Eigen::VectorXd weights =
    normalizeObservationWeights(observation_weights, y.rows());
  Eigen::MatrixXd theta(y.rows(), y.cols());

  for (Eigen::Index k = 0; k < y.cols(); ++k) {
    auto fit =
      uniformObservationWeights(weights)
        ? fitOls(x, Eigen::VectorXd(y.col(k)), intercept)
        : fitWeightedOls(x, Eigen::VectorXd(y.col(k)), intercept, weights);
    const auto& [beta0, beta] = fit;
    theta.col(k) = x * beta - y.col(k);
    if (intercept) {
      theta.col(k).array() += beta0;
      const double center = weightedMean(theta.col(k), weights)(0);
      theta.col(k).array() -= center;
    }
  }

  for (Eigen::Index i = 0; i < y.rows(); ++i) {
    if (weights(i) == 0.0) {
      theta.row(i).setZero();
    }
  }

  return theta;
}

template<typename T>
Eigen::MatrixXd
quadraticOlsDualPoint(
  const Eigen::MatrixBase<T>& x,
  const Eigen::MatrixXd& y,
  const Eigen::VectorXd& x_centers,
  const Eigen::VectorXd& x_scales,
  const JitNormalization jit_normalization,
  const bool intercept,
  const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
{
  const Eigen::VectorXd weights =
    normalizeObservationWeights(observation_weights, y.rows());
  Eigen::MatrixXd x_effective = x;
  const bool center = jit_normalization == JitNormalization::Center ||
                      jit_normalization == JitNormalization::Both;
  const bool scale = jit_normalization == JitNormalization::Scale ||
                     jit_normalization == JitNormalization::Both;

  for (Eigen::Index i = 0; i < weights.size(); ++i) {
    if (weights(i) == 0.0) {
      x_effective.row(i).setZero();
    }
  }
  if (center && !intercept) {
    x_effective.rowwise() -= x_centers.transpose();
  }
  // Mask again after centering so ignored rows cannot overflow during scaling.
  for (Eigen::Index i = 0; i < weights.size(); ++i) {
    if (weights(i) == 0.0) {
      x_effective.row(i).setZero();
    }
  }
  if (scale) {
    x_effective.array().rowwise() /= x_scales.transpose().array();
  }

  return quadraticOlsDualPointFromDesign(
    x_effective, y, intercept, observation_weights);
}

template<typename T>
Eigen::MatrixXd
quadraticOlsDualPoint(
  const Eigen::SparseMatrixBase<T>& x,
  const Eigen::MatrixXd& y,
  const Eigen::VectorXd& x_centers,
  const Eigen::VectorXd& x_scales,
  const JitNormalization jit_normalization,
  const bool intercept,
  const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
{
  const Eigen::VectorXd weights =
    normalizeObservationWeights(observation_weights, y.rows());
  const bool center = jit_normalization == JitNormalization::Center ||
                      jit_normalization == JitNormalization::Both;
  const bool scale = jit_normalization == JitNormalization::Scale ||
                     jit_normalization == JitNormalization::Both;

  Eigen::SparseMatrix<double> x_effective = x;
  for (int k = 0; k < x_effective.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(x_effective, k); it;
         ++it) {
      if (weights(it.row()) == 0.0) {
        it.valueRef() = 0.0;
      } else if (scale) {
        it.valueRef() /= x_scales(it.col());
      }
    }
  }

  if (!center || intercept) {
    return quadraticOlsDualPointFromDesign(
      x_effective, y, intercept, observation_weights);
  }

  Eigen::VectorXd effective_centers = x_centers;
  if (scale) {
    effective_centers.array() /= x_scales.array();
  }

  if (!uniformObservationWeights(observation_weights)) {
    // Without an intercept, explicit centering changes the constrained design.
    // Fit that design directly rather than adding a free intercept to sparse
    // QR.
    Eigen::MatrixXd centered =
      Eigen::MatrixXd(x_effective).rowwise() - effective_centers.transpose();
    for (Eigen::Index i = 0; i < weights.size(); ++i) {
      if (weights(i) == 0.0) {
        centered.row(i).setZero();
      }
    }
    return quadraticOlsDualPointFromDesign(
      centered, y, false, observation_weights);
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
computeDualFromPoint(
  const Eigen::VectorXd& beta,
  Eigen::MatrixXd theta,
  const std::unique_ptr<Loss>& loss,
  const SortedL1Norm& sl1_norm,
  const Eigen::ArrayXd& lambda,
  const MatrixType& x,
  const Eigen::MatrixXd& y,
  const Eigen::VectorXd& x_centers,
  const Eigen::VectorXd& x_scales,
  const JitNormalization& jit_normalization,
  const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
{
  const int n = x.rows();
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, n);
  Eigen::VectorXd gradient(beta.size());

  updateGradient(
    gradient, x, theta, x_centers, x_scales, weights, jit_normalization);

  const double dual_norm = sl1_norm.dualNorm(gradient, lambda);
  theta.array() /= std::max(1.0, dual_norm);

  return loss->dual(theta, y, weights);
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
            const bool intercept,
            const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
{
  Eigen::MatrixXd theta =
    loss->dualPoint(eta, y, intercept, observation_weights);
  return computeDualFromPoint(beta,
                              theta,
                              loss,
                              sl1_norm,
                              lambda,
                              x,
                              y,
                              x_centers,
                              x_scales,
                              jit_normalization,
                              observation_weights);
}

} // namespace slope
