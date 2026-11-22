#include <cassert>
#include <cmath>
#include <limits>
#include <slope/constants.h>
#include <slope/losses/poisson.h>
#include <slope/observation_weights.h>
#include <stdexcept>

namespace slope {

double
Poisson::loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  return (eta.array().exp() - y.array() * eta.array()).mean();
}

double
Poisson::loss(const Eigen::MatrixXd& eta,
              const Eigen::MatrixXd& y,
              const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  double value = 0.0;
  for (Eigen::Index i = 0; i < y.rows(); ++i) {
    if (weights(i) > 0.0) {
      for (Eigen::Index k = 0; k < y.cols(); ++k) {
        value += weights(i) * (std::exp(eta(i, k)) - y(i, k) * eta(i, k));
      }
    }
  }
  return value / y.size();
}

double
Poisson::dual(const Eigen::MatrixXd& theta,
              const Eigen::MatrixXd& y,
              const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  double value = 0.0;

  for (Eigen::Index i = 0; i < y.rows(); ++i) {
    if (weights(i) == 0.0) {
      continue;
    }
    for (Eigen::Index k = 0; k < y.cols(); ++k) {
      const double current = theta(i, k) + y(i, k);
      if (!std::isfinite(current) || current < 0.0) {
        throw std::domain_error(
          "Poisson dual means must be finite and nonnegative");
      }
      value += weights(i) *
               (current == 0.0 ? 0.0 : current * (1.0 - std::log(current)));
    }
  }

  return value / y.size();
}

Eigen::MatrixXd
Poisson::dualPoint(const Eigen::MatrixXd& eta,
                   const Eigen::MatrixXd& y,
                   const bool fit_intercept)
{
  return dualPoint(eta, y, fit_intercept, Eigen::VectorXd());
}

Eigen::MatrixXd
Poisson::dualPoint(const Eigen::MatrixXd& eta,
                   const Eigen::MatrixXd& y,
                   const bool fit_intercept,
                   const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  if (!fit_intercept) {
    return residualOnSupport(eta, y, weights);
  }
  if (eta.cols() != 1 || y.cols() != 1 || eta.rows() != y.rows()) {
    throw std::invalid_argument(
      "Poisson dual points require matching column vectors");
  }
  Eigen::MatrixXd safe_eta = eta;
  double maximum = -std::numeric_limits<double>::infinity();
  for (Eigen::Index i = 0; i < weights.size(); ++i) {
    if (weights(i) == 0.0) {
      safe_eta.row(i).setZero();
    } else {
      maximum = std::max(maximum, eta(i, 0));
    }
  }
  if (!safe_eta.allFinite()) {
    throw std::invalid_argument("Poisson linear predictors must be finite");
  }

  const double response_sum = weights.dot(y.col(0));
  if (response_sum <= 0.0) {
    return Eigen::MatrixXd::Zero(y.rows(), 1);
  }

  Eigen::VectorXd exponential = Eigen::VectorXd::Zero(y.rows());
  for (Eigen::Index i = 0; i < weights.size(); ++i) {
    if (weights(i) > 0.0) {
      exponential(i) = std::exp(eta(i, 0) - maximum);
    }
  }
  Eigen::MatrixXd mean =
    exponential * (response_sum / weights.dot(exponential));

  return domainSafePoint(
    mean, y, 0.0, std::numeric_limits<double>::infinity(), weights);
}

Eigen::MatrixXd
Poisson::hessianDiagonal(const Eigen::MatrixXd& eta)
{
  return eta.array().exp();
}

Eigen::MatrixXd
Poisson::preprocessResponse(const Eigen::MatrixXd& y)
{
  if ((y.array() < 0).any()) {
    throw std::invalid_argument("Response must be non-negative");
  }

  return y;
}

void
Poisson::updateIntercept(Eigen::VectorXd& beta0,
                         const Eigen::MatrixXd& eta,
                         const Eigen::MatrixXd& y)
{
  Eigen::VectorXd residual = this->residual(eta, y);
  double grad = residual.mean();
  double hess = eta.array().exp().mean();

  beta0(0) -= grad / hess;
}

Eigen::MatrixXd
Poisson::link(const Eigen::MatrixXd& mu)
{
  return mu.array().max(constants::P_MIN).log();
}

Eigen::MatrixXd
Poisson::inverseLink(const Eigen::MatrixXd& eta)
{
  return eta.array().exp();
}

Eigen::MatrixXd
Poisson::predict(const Eigen::MatrixXd& eta)
{
  return inverseLink(eta);
}

} // namespace slope
