#include <cassert>
#include <slope/math.h>
#include <slope/solvers/slope_threshold.h>

namespace slope {

std::tuple<double, int>
slopeThreshold(const double x,
               const int j,
               const Eigen::ArrayXd& lambda_cumsum,
               const Clusters& clusters)
{
  return slopeThreshold(x, 1.0, j, lambda_cumsum, clusters);
}

std::tuple<double, int>
slopeThreshold(const double gamma,
               const double omega,
               const int j,
               const Eigen::ArrayXd& lambda_cumsum,
               const Clusters& clusters)
{
  using std::size_t;

  assert(j >= 0 && j < clusters.size());
  assert(omega > 0);

  const auto& coeffs = clusters.coeffs();
  const auto& pointers = clusters.pointers();
  const int n_clusters = coeffs.size();
  const size_t ptr_j = pointers[j];
  const size_t cluster_size = pointers[j + 1] - pointers[j];
  const double abs_gamma = std::abs(gamma);
  const int sign_gamma = sign(gamma);

  // getLambdaSum(start, len) returns sum of lambdas from start to start+len-1
  auto getLambdaSum = [&](size_t start, size_t len) -> double {
    return lambda_cumsum(start + len) - lambda_cumsum(start);
  };

  // Determine whether the update moves upward.
  const double current_lambda_sum = getLambdaSum(ptr_j, cluster_size);
  const bool direction_up = abs_gamma - current_lambda_sum > omega * coeffs[j];

  if (direction_up) {
    double lo = current_lambda_sum;

    for (int k = j - 1; k >= 0; --k) {
      const double c_k = coeffs[k];

      if (abs_gamma - lo < omega * c_k) {
        return { (gamma - sign_gamma * lo) / omega, k + 1 };
      }

      const size_t start = pointers[k];
      const double hi = getLambdaSum(start, cluster_size);

      if (abs_gamma - hi <= omega * c_k) {
        return { sign_gamma * c_k, k };
      }

      lo = hi;
    }

    return { (gamma - sign_gamma * lo) / omega, 0 };
  } else {
    // Moving down in the cluster ordering
    double hi = current_lambda_sum;

    for (int k = j + 1; k < n_clusters; ++k) {
      const size_t end = pointers[k + 1];

      const double c_k = coeffs[k];

      if (abs_gamma > hi + omega * c_k) {
        return { (gamma - sign_gamma * hi) / omega, k - 1 };
      }

      const double lo = getLambdaSum(end - cluster_size, cluster_size);

      if (abs_gamma >= lo + omega * c_k) {
        return { sign_gamma * c_k, k };
      }

      hi = lo;
    }

    if (abs_gamma > hi) {
      return { (gamma - sign_gamma * hi) / omega, n_clusters - 1 };
    } else {
      // Zero cluster case
      return { 0, n_clusters };
    }
  }
}

}
