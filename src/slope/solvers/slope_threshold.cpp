#include "slope_threshold.h"
#include "../math.h"
#include <algorithm>

namespace slope {

std::tuple<double, int>
slopeThreshold(const double x,
               const int j,
               const Eigen::ArrayXd lambdas,
               const Clusters& clusters)
{
  using std::size_t;

  const size_t cluster_size = clusters.cluster_size(j);
  const double abs_x = std::abs(x);
  const int sign_x = sign(x);
  const size_t n_lambda = lambdas.size();

  // Prepare a lazy cumulative sum of lambdas.
  // cum[i] holds sum_{k=0}^{i-1} lambdas(k) with cum[0] = 0.
  std::vector<double> cum(n_lambda + 1, 0.0);
  size_t computed = 0; // Last index for which cum has been computed.

  // getCum(i) computes and returns cum[i] on demand.
  auto getCum = [&](size_t i) -> double {
    i = std::min(i, n_lambda); // Prevent out-of-bounds access
    while (computed < i) {
      computed++;
      cum[computed] = cum[computed - 1] + lambdas(computed - 1);
    }
    return cum[i];
  };

  // getLambdaSum(start, len) returns sum of lambdas from start to start+len-1
  auto getLambdaSum = [&](size_t start, size_t len) -> double {
    return getCum(start + len) - getCum(start);
  };

  // Determine whether the update moves upward.
  int ptr_j = clusters.pointer(j);
  size_t len_j = std::min(n_lambda - ptr_j, cluster_size);
  const bool direction_up =
    abs_x - getLambdaSum(ptr_j, len_j) > clusters.coeff(j);

  if (direction_up) {
    // Moving up in the cluster ordering
    size_t start = clusters.pointer(j + 1);
    size_t len = std::min(n_lambda - start, cluster_size);
    double lo = start < n_lambda ? getLambdaSum(start, len) : 0.0;

    for (int k = j; k >= 0; --k) {
      start = clusters.pointer(k);
      len = std::min(n_lambda - start, cluster_size);
      double hi = getLambdaSum(start, len);
      double c_k = clusters.coeff(k);

      if (abs_x < lo + c_k) {
        return { x - sign_x * lo, k + 1 };
      } else if (abs_x <= hi + c_k) {
        return { sign_x * c_k, k };
      }
      lo = hi;
    }

    return { x - sign_x * lo, 0 };
  } else {
    // Moving down in the cluster ordering
    int end = clusters.pointer(j + 1);
    double hi = getLambdaSum(end - cluster_size, cluster_size);

    for (int k = j + 1; k < clusters.n_clusters(); ++k) {
      end = clusters.pointer(k + 1);
      double lo = getLambdaSum(end - cluster_size, cluster_size);
      double c_k = clusters.coeff(k);

      if (abs_x > hi + c_k) {
        return { x - sign_x * hi, k - 1 };
      } else if (abs_x >= lo + c_k) {
        return { sign_x * c_k, k };
      }
      hi = lo;
    }

    if (abs_x > hi) {
      return { x - sign_x * hi, clusters.n_clusters() - 1 };
    } else {
      return { 0, clusters.n_clusters() - 1 };
    }
  }
}

}
