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

bool
checkKktViolations(Eigen::VectorXd& gradient,
                   const Eigen::VectorXd& beta,
                   const Eigen::ArrayXd& lambda_curr,
                   std::vector<int>& working_set,
                   const std::vector<int>& strong_set,
                   const std::vector<int>& full_set,
                   const Eigen::MatrixXd& x,
                   const Eigen::MatrixXd& residual,
                   const Eigen::VectorXd& x_centers,
                   const Eigen::VectorXd& x_scales,
                   JitNormalization jit_normalization)
{
  // Check for violations in the strong set first
  updateGradient(gradient,
                 x,
                 residual,
                 strong_set,
                 x_centers,
                 x_scales,
                 Eigen::VectorXd::Ones(x.rows()),
                 jit_normalization);

  auto violations =
    setDiff(kktCheck(gradient, beta, lambda_curr, strong_set), working_set);

  if (violations.empty()) {
    // Now check for violations in the full set
    updateGradient(gradient,
                   x,
                   residual,
                   full_set,
                   x_centers,
                   x_scales,
                   Eigen::VectorXd::Ones(x.rows()),
                   jit_normalization);

    violations =
      setDiff(kktCheck(gradient, beta, lambda_curr, full_set), working_set);

    if (violations.empty()) {
      return true; // No violations found
    } else {
      working_set = setUnion(working_set, violations);
      return false; // Violations found and working set updated
    }
  } else {
    working_set = setUnion(working_set, violations);
    return false; // Violations found and working set updated
  }
}

} // namespace slope
