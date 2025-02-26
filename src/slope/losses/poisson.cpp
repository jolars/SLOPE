#include "poisson.h"
#include "../math.h"

namespace slope {

double
Poisson::loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  return -(y.array() * eta.array() - eta.array().exp() - logFactorial(y))
            .mean();
}

double
Poisson::dual(const Eigen::MatrixXd& theta,
              const Eigen::MatrixXd& y,
              const Eigen::VectorXd&)
{
  const Eigen::ArrayXd e = theta + y;

  assert(theta.allFinite() && "theta is not finite");
  assert((e >= 0).all() &&
         "Dual function is not defined for negative residuals");

  return (e * (1.0 - e.log()) + logFactorial(y)).mean();
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
  return mu.array().log();
}

} // namespace slope
