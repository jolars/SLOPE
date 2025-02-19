#include "quadratic.h"

namespace slope {

double
Quadratic::loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  return (eta - y).squaredNorm() / (2.0 * y.rows());
}

double
Quadratic::dual(const Eigen::MatrixXd& theta,
                const Eigen::MatrixXd& y,
                const Eigen::VectorXd&)
{
  const int n = y.rows();

  return (y.squaredNorm() - (theta + y).squaredNorm()) / (2.0 * n);
}

Eigen::MatrixXd
Quadratic::residual(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  return eta - y;
}

Eigen::MatrixXd
Quadratic::preprocessResponse(const Eigen::MatrixXd& y)
{
  return y;
}

void
Quadratic::updateWeightsAndWorkingResponse(Eigen::VectorXd& w,
                                           Eigen::VectorXd& z,
                                           const Eigen::VectorXd&,
                                           const Eigen::VectorXd& y)
{
  w.setOnes();
  z = y;
}

Eigen::MatrixXd
Quadratic::link(const Eigen::MatrixXd& mu)
{
  return mu;
}

// double
// Quadratic::nullDeviance(const Eigen::MatrixXd& y, const bool intercept)
// {
//   double beta0 = intercept ? y.mean() : 0.0;
//
//   Eigen::MatrixXd eta(y.rows(), y.cols());
//   eta.setConstant(beta0);
//
//   return deviance(eta, y);
// }

} // namespace slope
