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
Quadratic::hessianDiagonal(const Eigen::MatrixXd& eta)
{
  return Eigen::MatrixXd::Ones(eta.rows(), eta.cols());
}

Eigen::MatrixXd
Quadratic::preprocessResponse(const Eigen::MatrixXd& y)
{
  return y;
}

void
Quadratic::updateWeightsAndWorkingResponse(Eigen::MatrixXd& w,
                                           Eigen::MatrixXd& z,
                                           const Eigen::MatrixXd&,
                                           const Eigen::MatrixXd& y)
{
  // Do nothing since weights are already one and working response is y.
}

Eigen::MatrixXd
Quadratic::link(const Eigen::MatrixXd& mu)
{
  return mu;
}

Eigen::MatrixXd
Quadratic::inverseLink(const Eigen::MatrixXd& eta)
{
  return eta;
}

Eigen::MatrixXd
Quadratic::predict(const Eigen::MatrixXd& eta)
{
  return eta;
}

} // namespace slope
