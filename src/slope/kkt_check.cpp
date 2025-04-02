#include "kkt_check.h"
#include "math.h"
#include "utils.h"
#include <Eigen/Core>

namespace slope {

std::vector<int>
kktCheck(const Eigen::VectorXd& gradient,
         const Eigen::VectorXd& beta,
         const Eigen::ArrayXd& lambda,
         const std::vector<int>& indices)
{
  using namespace Eigen;

  std::vector<int> out;

  int pm = beta.size();

  if (pm == 0) {
    return out;
  }

  VectorXd flat_abs_gradient = gradient(indices).cwiseAbs();

  auto ord = sortIndex(flat_abs_gradient, true);

  // TODO: Use dualNorm() instead
  ArrayXd diff = flat_abs_gradient(ord).array() - lambda.head(indices.size());
  ArrayXb tmp = cumSum(diff) >= 0.0;
  inversePermute(tmp, ord);

  auto which_violations = which(tmp);

  for (const auto& ind : which_violations) {
    out.emplace_back(indices[ind]);
  }

  return out;
}

} // namespace slope
