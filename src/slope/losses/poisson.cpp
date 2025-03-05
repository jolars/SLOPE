#include "poisson.h"
#include "slope/constants.h"

namespace slope {

double
Poisson::loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  return -(y.array() * eta.array() - eta.array().exp()).mean();
}

double
Poisson::dual(const Eigen::MatrixXd& theta,
              const Eigen::MatrixXd& y,
              const Eigen::VectorXd&)
{
  const Eigen::ArrayXd e = theta + y;

  assert(theta.allFinite() && "theta is not finite");

  return (e * (1.0 - e.max(constants::P_MIN).log())).mean();
}

Eigen::MatrixXd
Poisson::residual(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  return eta.array().exp() - y.array();
}

void
Poisson::updateWeightsAndWorkingResponse(Eigen::VectorXd& w,
                                         Eigen::VectorXd& z,
                                         const Eigen::VectorXd& eta,
                                         const Eigen::VectorXd& y)
{
  w = eta.array().exp();
  z = eta.array() - 1.0 + y.array() / w.array();
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
