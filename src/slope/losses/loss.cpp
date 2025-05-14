#include "loss.h"
#include <Eigen/Core>

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

} // namespace slope
