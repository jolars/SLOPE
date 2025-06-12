#include "logistic.h"
#include "../constants.h"
#include "../math.h"

namespace slope {

double
Logistic::loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  double loss =
    eta.array().exp().log1p().sum() - y.reshaped().dot(eta.reshaped());
  return loss / y.rows();
}

double
Logistic::dual(const Eigen::MatrixXd& theta,
               const Eigen::MatrixXd& y,
               const Eigen::VectorXd&)
{
  int n = y.rows();

  Eigen::VectorXd eta = link(theta + y);

  return loss(eta, y) - theta.reshaped().dot(eta) / n;
}

Eigen::MatrixXd
Logistic::hessianDiagonal(const Eigen::MatrixXd& eta)
{
  const auto pr = inverseLink(eta);
  return pr.array() * (1.0 - pr.array());
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
  Eigen::ArrayXd pr =
    1.0 / (1.0 + (-eta).array().min(constants::MAX_EXP).exp());

  return pr.max(constants::P_MIN).min(constants::P_MAX);
}

Eigen::MatrixXd
Logistic::predict(const Eigen::MatrixXd& eta)
{
  Eigen::MatrixXd prob = inverseLink(eta);
  return prob.unaryExpr([](double pr) { return pr > 0.5 ? 1.0 : 0.0; });
}

} // namespace slope
