#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <limits>
#include <slope/losses/loss.h>
#include <slope/observation_weights.h>
#include <stdexcept>

namespace slope {

double
Loss::loss(const Eigen::MatrixXd& eta,
           const Eigen::MatrixXd& y,
           const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  if (!detail::uniformObservationWeights(weights)) {
    throw std::invalid_argument(
      "This loss does not support observation weights");
  }
  return loss(eta, y);
}

void
Loss::updateWeightsAndWorkingResponse(Eigen::MatrixXd& w,
                                      Eigen::MatrixXd& z,
                                      const Eigen::MatrixXd& eta,
                                      const Eigen::MatrixXd& y)
{
  w = hessianDiagonal(eta);
  z = eta.array() + (y.array() - inverseLink(eta).array()) / w.array();
}

void
Loss::updateWeightsAndWorkingResponse(
  Eigen::MatrixXd& w,
  Eigen::MatrixXd& z,
  const Eigen::MatrixXd& eta,
  const Eigen::MatrixXd& y,
  const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  Eigen::MatrixXd safe_eta = eta;
  Eigen::MatrixXd safe_y = y;
  for (Eigen::Index i = 0; i < weights.size(); ++i) {
    if (weights(i) == 0.0) {
      safe_eta.row(i).setZero();
      safe_y.row(i).setZero();
    }
  }

  w = Eigen::MatrixXd::Ones(eta.rows(), eta.cols());
  z = safe_y;
  updateWeightsAndWorkingResponse(w, z, safe_eta, safe_y);
  for (Eigen::Index i = 0; i < weights.size(); ++i) {
    if (weights(i) == 0.0) {
      w.row(i).setZero();
      z.row(i) = eta.row(i);
    } else {
      w.row(i) *= weights(i);
    }
  }
}

double
Loss::deviance(const Eigen::MatrixXd& eta,
               const Eigen::MatrixXd& y,
               const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  if (detail::uniformObservationWeights(weights)) {
    return deviance(eta, y);
  }
  return 2.0 * (loss(eta, y, weights) - loss(link(y), y, weights));
}

Eigen::MatrixXd
Loss::residual(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  return inverseLink(eta) - y;
}

Eigen::MatrixXd
Loss::dualPoint(const Eigen::MatrixXd& eta,
                const Eigen::MatrixXd& y,
                const bool fit_intercept)
{
  Eigen::MatrixXd theta = residual(eta, y);
  if (fit_intercept) {
    theta.rowwise() -= theta.colwise().mean();
  }
  return theta;
}

Eigen::MatrixXd
Loss::dualPoint(const Eigen::MatrixXd& eta,
                const Eigen::MatrixXd& y,
                const bool fit_intercept,
                const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  if (detail::uniformObservationWeights(weights)) {
    return dualPoint(eta, y, fit_intercept);
  }
  Eigen::MatrixXd theta = residualOnSupport(eta, y, weights);
  if (fit_intercept) {
    const Eigen::RowVectorXd center = detail::weightedMean(theta, weights);
    for (Eigen::Index i = 0; i < weights.size(); ++i) {
      if (weights(i) > 0.0) {
        theta.row(i) -= center;
      }
    }
  }
  return theta;
}

Eigen::MatrixXd
Loss::residualOnSupport(const Eigen::MatrixXd& eta,
                        const Eigen::MatrixXd& y,
                        const Eigen::VectorXd& weights)
{
  Eigen::MatrixXd safe_eta = eta;
  for (Eigen::Index i = 0; i < weights.size(); ++i) {
    if (weights(i) == 0.0) {
      safe_eta.row(i).setZero();
    }
  }
  Eigen::MatrixXd theta = residual(safe_eta, y);
  for (Eigen::Index i = 0; i < weights.size(); ++i) {
    if (weights(i) == 0.0) {
      theta.row(i).setZero();
    }
  }
  return theta;
}

Eigen::MatrixXd
Loss::domainSafePoint(const Eigen::MatrixXd& mean,
                      const Eigen::MatrixXd& y,
                      const double lower,
                      const double upper) const
{
  return domainSafePoint(mean, y, lower, upper, Eigen::VectorXd());
}

Eigen::MatrixXd
Loss::domainSafePoint(const Eigen::MatrixXd& mean,
                      const Eigen::MatrixXd& y,
                      const double lower,
                      const double upper,
                      const Eigen::VectorXd& observation_weights) const
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  Eigen::MatrixXd centered = mean - y;
  centered.rowwise() -= detail::weightedMean(centered, weights);

  const double response_mean = detail::weightedMean(y, weights)(0);
  double step = 1.0;

  for (Eigen::Index i = 0; i < centered.size(); ++i) {
    if (weights(i % y.rows()) == 0.0) {
      continue;
    }
    const double target = y(i) + centered(i);
    const double direction = target - response_mean;
    if (direction < 0.0) {
      step = std::min(step, (response_mean - lower) / -direction);
    } else if (direction > 0.0 && std::isfinite(upper)) {
      step = std::min(step, (upper - response_mean) / direction);
    }
  }

  step = std::clamp(step, 0.0, 1.0);
  if (step > 0.0 && step < 1.0) {
    step *= 1.0 - std::sqrt(std::numeric_limits<double>::epsilon());
  }

  Eigen::MatrixXd theta(centered.rows(), centered.cols());
  for (Eigen::Index i = 0; i < centered.size(); ++i) {
    if (weights(i % y.rows()) == 0.0) {
      theta(i) = 0.0;
      continue;
    }
    const double anchor = response_mean - y(i);
    theta(i) = anchor + step * (centered(i) - anchor);
  }
  return theta;
}

} // namespace slope
