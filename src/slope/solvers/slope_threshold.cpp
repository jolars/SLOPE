#include "slope_threshold.h"
#include "../math.h"

namespace slope {

std::tuple<double, int>
slopeThreshold(const double x,
               const int j,
               const Eigen::ArrayXd& lambda_cumsum,
               const Clusters& clusters)
{
  using std::size_t;

  assert(j >= 0 && j < clusters.size());

  const size_t cluster_size = clusters.cluster_size(j);
  const double abs_x = std::abs(x);
  const int sign_x = sign(x);

  // getLambdaSum(start, len) returns sum of lambdas from start to start+len-1
  auto getLambdaSum = [&](size_t start, size_t len) -> double {
    return lambda_cumsum(start + len) - lambda_cumsum(start);
  };

  // Determine whether the update moves upward.
  int ptr_j = clusters.pointer(j);
  const bool direction_up =
    abs_x - getLambdaSum(ptr_j, cluster_size) > clusters.coeff(j);

  if (direction_up) {
    size_t start = clusters.pointer(j);
    double lo = getLambdaSum(start, cluster_size);

    for (int k = j - 1; k >= 0; --k) {
      double c_k = clusters.coeff(k);

      if (abs_x - lo < c_k && k < j) {
        return { x - sign_x * lo, k + 1 };
      }

      start = clusters.pointer(k);
      double hi = getLambdaSum(start, cluster_size);

      if (abs_x - hi <= c_k) {
        return { sign_x * c_k, k };
      }

      lo = hi;
    }

    return { x - sign_x * lo, 0 };
  } else {
    // Moving down in the cluster ordering
    int end = clusters.pointer(j + 1);
    double hi = getLambdaSum(end - cluster_size, cluster_size);

    for (int k = j + 1; k < clusters.size(); ++k) {
      end = clusters.pointer(k + 1);

      double c_k = clusters.coeff(k);

      if (abs_x > hi + c_k) {
        return { x - sign_x * hi, k - 1 };
      }

      double lo = getLambdaSum(end - cluster_size, cluster_size);

      if (abs_x >= lo + c_k) {
        return { sign_x * c_k, k };
      }

      hi = lo;
    }

    if (abs_x > hi) {
      return { x - sign_x * hi, clusters.size() - 1 };
    } else {
      // Zero cluster case
      return { 0, clusters.size() };
    }
  }
}

}
