#include "slope_threshold.h"
#include "../math.h"

namespace slope {

std::tuple<double, int>
slopeThreshold(const double x,
               const int j,
               const Eigen::ArrayXd lambdas,
               const Clusters& clusters)
{
  const Eigen::Index cluster_size = clusters.cluster_size(j);
  const double abs_x = std::abs(x);
  const int sign_x = sign(x);

  const bool direction_up =
    abs_x - lambdas.segment(clusters.pointer(j), cluster_size).sum() >
    clusters.coeff(j);

  // TODO: These partial lambda sums may actually overlap, which
  // means that we are doing a lot of redundant computations in the case when
  // the cluster is large. We could do better here.

  if (direction_up) {
    int start = clusters.pointer(j + 1);
    int end = std::min(lambdas.size() - start, cluster_size);
    double lo =
      start < lambdas.size() ? lambdas.segment(start, end).sum() : 0.0;

    for (int k = j; k >= 0; k--) {
      start = clusters.pointer(k);
      Eigen::Index end = std::min(lambdas.size() - start, cluster_size);
      double hi = lambdas.segment(start, end).sum();
      double c_k = clusters.coeff(k);

      if (abs_x < lo + c_k) {
        // We are in-between clusters.
        return { x - sign_x * lo, k + 1 };
      } else if (abs_x <= hi + c_k) {
        // We are in a cluster.
        return { sign_x * c_k, k };
      }
      lo = hi;
    }

    return { x - sign_x * lo, 0 };
  } else {
    int end = clusters.pointer(j + 1);
    double hi = lambdas.segment(end - cluster_size, cluster_size).sum();

    for (int k = j + 1; k < clusters.n_clusters(); ++k) {
      end = clusters.pointer(k + 1);
      double lo = lambdas.segment(end - cluster_size, cluster_size).sum();
      double c_k = clusters.coeff(k);

      if (abs_x > hi + c_k) {
        // We are in-between clusters.
        return { x - sign_x * hi, k - 1 };
      } else if (abs_x >= lo + c_k) {
        // We are in a cluster.
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
