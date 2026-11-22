#include <cmath>
#include <limits>
#include <slope/constants.h>
#include <slope/losses/logistic.h>
#include <slope/math.h>
#include <slope/observation_weights.h>
#include <stdexcept>

namespace {

double
softplus(const double x)
{
  return x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

double
stableSigmoid(const double x)
{
  if (x >= 0.0) {
    return 1.0 / (1.0 + std::exp(-x));
  }
  const double exp_x = std::exp(x);
  return exp_x / (1.0 + exp_x);
}

double
xLogX(const double x)
{
  return x == 0.0 ? 0.0 : x * std::log(x);
}

} // namespace

namespace slope {

double
Logistic::loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  const double value =
    eta.unaryExpr(&softplus).sum() - y.reshaped().dot(eta.reshaped());
  return value / y.rows();
}

double
Logistic::loss(const Eigen::MatrixXd& eta,
               const Eigen::MatrixXd& y,
               const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  double value = 0.0;
  for (Eigen::Index i = 0; i < y.rows(); ++i) {
    if (weights(i) > 0.0) {
      for (Eigen::Index k = 0; k < y.cols(); ++k) {
        value += weights(i) * (softplus(eta(i, k)) - y(i, k) * eta(i, k));
      }
    }
  }
  return value / y.rows();
}

double
Logistic::dual(const Eigen::MatrixXd& theta,
               const Eigen::MatrixXd& y,
               const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  double entropy = 0.0;

  for (Eigen::Index i = 0; i < y.rows(); ++i) {
    if (weights(i) == 0.0) {
      continue;
    }
    for (Eigen::Index k = 0; k < y.cols(); ++k) {
      const double value = theta(i, k) + y(i, k);
      if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
        throw std::domain_error(
          "Logistic dual means must be finite and lie in [0, 1]");
      }
      entropy += weights(i) * (xLogX(value) + xLogX(1.0 - value));
    }
  }

  return -entropy / y.rows();
}

Eigen::MatrixXd
Logistic::dualPoint(const Eigen::MatrixXd& eta,
                    const Eigen::MatrixXd& y,
                    const bool fit_intercept)
{
  return dualPoint(eta, y, fit_intercept, Eigen::VectorXd());
}

Eigen::MatrixXd
Logistic::dualPoint(const Eigen::MatrixXd& eta,
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
      "Logistic dual points require matching column vectors");
  }
  Eigen::MatrixXd safe_eta = eta;
  for (Eigen::Index i = 0; i < weights.size(); ++i) {
    if (weights(i) == 0.0) {
      safe_eta.row(i).setZero();
    }
  }
  if (!safe_eta.allFinite()) {
    throw std::invalid_argument("Logistic linear predictors must be finite");
  }

  const double response_mean = detail::weightedMean(y, weights)(0);
  if (response_mean <= 0.0 || response_mean >= 1.0) {
    return domainSafePoint(
      Eigen::MatrixXd::Constant(y.rows(), 1, response_mean),
      y,
      0.0,
      1.0,
      weights);
  }

  const double target = std::log(response_mean) - std::log1p(-response_mean);
  double lower = target - safe_eta.maxCoeff();
  double upper = target - safe_eta.minCoeff();
  double shift = std::clamp(0.0, lower, upper);
  Eigen::MatrixXd mean(eta.rows(), 1);

  for (int iteration = 0; iteration < 40; ++iteration) {
    mean = inverseLink((safe_eta.array() + shift).matrix());
    const double gradient = weights.dot((mean - y).col(0));
    if (std::abs(gradient) <=
        16.0 * std::numeric_limits<double>::epsilon() * y.rows()) {
      break;
    }

    if (gradient < 0.0) {
      lower = shift;
    } else {
      upper = shift;
    }

    const double hessian =
      (weights.array() * mean.col(0).array() * (1.0 - mean.col(0).array()))
        .sum();
    const double newton = shift - gradient / hessian;
    shift = std::isfinite(newton) && newton > lower && newton < upper
              ? newton
              : 0.5 * (lower + upper);
  }

  mean = inverseLink((safe_eta.array() + shift).matrix());
  return domainSafePoint(mean, y, 0.0, 1.0, weights);
}

Eigen::MatrixXd
Logistic::hessianDiagonal(const Eigen::MatrixXd& eta)
{
  const auto pr = inverseLink(eta);
  constexpr double min_hessian = constants::P_MIN * (1.0 - constants::P_MIN);
  return (pr.array() * (1.0 - pr.array())).max(min_hessian);
}

Eigen::MatrixXd
Logistic::preprocessResponse(const Eigen::MatrixXd& y)
{
  // Check if the response is in {0, 1} and convert it otherwise
  Eigen::MatrixXd y_clamped = y.array().min(1.0).max(0.0);

  // Throw an error if the response is not binary
  if ((y_clamped.array() != 0.0 && y_clamped.array() != 1.0).any()) {
    throw std::invalid_argument("Response must be binary");
  }

  return y_clamped;
}

Eigen::MatrixXd
Logistic::link(const Eigen::MatrixXd& mu)
{
  return mu.unaryExpr([](const double& x) {
    return logit(std::clamp(x, constants::P_MIN, constants::P_MAX));
  });
}

Eigen::MatrixXd
Logistic::inverseLink(const Eigen::MatrixXd& eta)
{
  return eta.unaryExpr([](const double value) { return stableSigmoid(value); });
}

Eigen::MatrixXd
Logistic::predict(const Eigen::MatrixXd& eta)
{
  Eigen::MatrixXd prob = inverseLink(eta);
  return prob.unaryExpr([](double pr) { return pr > 0.5 ? 1.0 : 0.0; });
}

} // namespace slope
