#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <limits>
#include <slope/losses/loss.h>

namespace slope {

void
Loss::updateWeightsAndWorkingResponse(Eigen::MatrixXd& w,
                                      Eigen::MatrixXd& z,
                                      const Eigen::MatrixXd& eta,
                                      const Eigen::MatrixXd& y)
{
  w = hessianDiagonal(eta);
  z = eta.array() + (y.array() - inverseLink(eta).array()) / w.array();
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
Loss::domainSafePoint(const Eigen::MatrixXd& mean,
                      const Eigen::MatrixXd& y,
                      const double lower,
                      const double upper) const
{
  Eigen::MatrixXd centered = mean - y;
  centered.rowwise() -= centered.colwise().mean();

  const double response_mean = y.mean();
  double step = 1.0;

  for (Eigen::Index i = 0; i < centered.size(); ++i) {
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
    const double anchor = response_mean - y(i);
    theta(i) = anchor + step * (centered(i) - anchor);
  }
  return theta;
}

} // namespace slope
