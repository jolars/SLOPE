#include "kkt_check.h"
#include "slope/math.h"
#include "slope/utils.h"
#include <Eigen/Core>

namespace slope {

std::vector<int>
kktCheck(const Eigen::MatrixXd& gradient,
         const Eigen::MatrixXd& beta,
         const Eigen::ArrayXd& lambda,
         const std::vector<int>& indices)
{
  using namespace Eigen;

  std::vector<int> out;

  int p = beta.rows();
  int m = beta.cols();

  if (p == 0) {
    return out;
  }

  VectorXd flat_abs_gradient = gradient(indices, all).reshaped().cwiseAbs();

  auto ord = sortIndex(flat_abs_gradient, true);

  // TODO: Use dualNorm() instead
  ArrayXd diff =
    flat_abs_gradient(ord).array() - lambda.head(flat_abs_gradient.size());
  ArrayXXb tmp = cumSum(diff).array() >= 0.0;
  tmp(ord, 0) = tmp.col(0).eval();
  tmp.resize(indices.size(), m);

  ArrayXb violations = tmp.rowwise().any();
  auto which_violations = which(violations);

  for (const auto& ind : which_violations) {
    out.emplace_back(indices[ind]);
  }

  return out;
}

}
