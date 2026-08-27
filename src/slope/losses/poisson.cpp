#include <cassert>
#include <cmath>
#include <limits>
#include <slope/constants.h>
#include <slope/losses/poisson.h>
#include <stdexcept>

namespace slope {

double
Poisson::loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  return (eta.array().exp() - y.array() * eta.array()).mean();
}

double
Poisson::dual(const Eigen::MatrixXd& theta,
              const Eigen::MatrixXd& y,
              const Eigen::VectorXd&)
{
  const Eigen::ArrayXXd mean = theta.array() + y.array();
  double value = 0.0;

  for (Eigen::Index i = 0; i < mean.size(); ++i) {
    const double current = mean(i);
    if (!std::isfinite(current) || current < 0.0) {
      throw std::domain_error(
        "Poisson dual means must be finite and nonnegative");
    }
    value += current == 0.0 ? 0.0 : current * (1.0 - std::log(current));
  }

  return value / y.rows();
}

Eigen::MatrixXd
Poisson::dualPoint(const Eigen::MatrixXd& eta,
                   const Eigen::MatrixXd& y,
                   const bool fit_intercept)
{
  if (!fit_intercept) {
    return residual(eta, y);
  }
  if (eta.cols() != 1 || y.cols() != 1 || eta.rows() != y.rows()) {
    throw std::invalid_argument(
      "Poisson dual points require matching column vectors");
  }
  if (!eta.allFinite()) {
    throw std::invalid_argument("Poisson linear predictors must be finite");
  }

  const double response_sum = y.sum();
  if (response_sum <= 0.0) {
    return Eigen::MatrixXd::Zero(y.rows(), 1);
  }

  const double maximum = eta.maxCoeff();
  Eigen::ArrayXd weights = (eta.array() - maximum).exp();
  Eigen::MatrixXd mean = (weights * (response_sum / weights.sum())).matrix();

  return domainSafePoint(mean, y, 0.0, std::numeric_limits<double>::infinity());
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
