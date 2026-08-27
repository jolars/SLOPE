#include "kkt_check.h"
#include <Eigen/Core>
#include <slope/math.h>
#include <slope/utils.h>

namespace slope {

namespace {

template<typename IndexAt>
std::vector<int>
kktCheckImpl(const Eigen::VectorXd& gradient,
             const Eigen::VectorXd& beta,
             const Eigen::ArrayXd& lambda,
             const int index_count,
             const IndexAt& index_at)
{
  using namespace Eigen;

  std::vector<int> out;

  if (beta.size() == 0) {
    return out;
  }

  ArrayXd abs_gradient(index_count);
  for (int i = 0; i < index_count; ++i) {
    abs_gradient(i) = std::abs(gradient(index_at(i)));
  }

  auto ord = sortIndex(abs_gradient, true);
  permute(abs_gradient, ord);

  ArrayXd diff = abs_gradient - lambda.head(index_count);
  ArrayXb tmp = cumSum(diff) >= 0.0;

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
    out.emplace_back(index_at(ord[i]));
  }

  std::sort(out.begin(), out.end());

  return out;
}

} // namespace

std::vector<int>
kktCheck(const Eigen::VectorXd& gradient,
         const Eigen::VectorXd& beta,
         const Eigen::ArrayXd& lambda,
         const std::vector<int>& indices)
{
  return kktCheckImpl(gradient,
                      beta,
                      lambda,
                      static_cast<int>(indices.size()),
                      [&indices](int i) { return indices[i]; });
}

std::vector<int>
kktCheck(const Eigen::VectorXd& gradient,
         const Eigen::VectorXd& beta,
         const Eigen::ArrayXd& lambda)
{
  return kktCheckImpl(
    gradient, beta, lambda, static_cast<int>(gradient.size()), [](int i) {
      return i;
    });
}

} // namespace slope
