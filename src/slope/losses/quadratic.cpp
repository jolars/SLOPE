#include <slope/losses/quadratic.h>
#include <slope/observation_weights.h>

namespace slope {

double
Quadratic::loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  return (eta - y).squaredNorm() / (2.0 * y.rows());
}

double
Quadratic::loss(const Eigen::MatrixXd& eta,
                const Eigen::MatrixXd& y,
                const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  double value = 0.0;
  for (Eigen::Index i = 0; i < y.rows(); ++i) {
    if (weights(i) > 0.0) {
      value += weights(i) * (eta.row(i) - y.row(i)).squaredNorm();
    }
  }
  return value / (2.0 * y.rows());
}

double
Quadratic::dual(const Eigen::MatrixXd& theta,
                const Eigen::MatrixXd& y,
                const Eigen::VectorXd& observation_weights)
{
  const Eigen::VectorXd weights =
    detail::normalizeObservationWeights(observation_weights, y.rows());
  double value = 0.0;
  for (Eigen::Index i = 0; i < y.rows(); ++i) {
    if (weights(i) > 0.0) {
      value -= weights(i) *
               (theta.row(i).squaredNorm() + 2.0 * theta.row(i).dot(y.row(i)));
    }
  }
  return value / (2.0 * y.rows());
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
  w = Eigen::MatrixXd::Ones(y.rows(), y.cols());
  z = y;
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
