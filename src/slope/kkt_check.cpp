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

  ArrayXd abs_gradient = gradient(indices).cwiseAbs();
  auto ord = sortIndex(abs_gradient, true);
  permute(abs_gradient, ord);

  ArrayXd diff = abs_gradient - lambda.head(indices.size());
  ArrayXb tmp = cumSum(diff) >= 0.0;

  // Find the last position where cumulative sum is non-negative
  int k = 0;
  if (tmp.size() > 0) {
    for (int i = tmp.size() - 1; i >= 0; --i) {
      if (tmp[i]) {
        k = i + 1;
        break;
      }
    }
  }

  out.reserve(k);
  for (int i = 0; i < k; ++i) {
    out.emplace_back(indices[ord[i]]);
  }

  std::sort(out.begin(), out.end());

  return out;
}

} // namegspace slope
