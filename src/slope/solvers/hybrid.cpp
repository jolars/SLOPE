/**
 * @file
 * @brief Hybrid solver implementation for SLOPE
 */

#include "hybrid.h"
#include "../sorted_l1_norm.h"
#include "math.h"
#include "slope/clusters.h"
#include "slope/losses/loss.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <memory>

namespace slope {
namespace solvers {

// Override for dense matrices
void
Hybrid::run(Eigen::VectorXd& beta0,
            Eigen::MatrixXd& beta,
            Eigen::MatrixXd& eta,
            Clusters& clusters,
            const Eigen::ArrayXd& lambda,
            const std::unique_ptr<Loss>& loss,
            const SortedL1Norm& penalty,
            const Eigen::MatrixXd& gradient,
            const std::vector<int>& working_set,
            const Eigen::MatrixXd& x,
            const Eigen::VectorXd& x_centers,
            const Eigen::VectorXd& x_scales,
            const Eigen::MatrixXd& y)
{
  runImpl(beta0,
          beta,
          eta,
          clusters,
          lambda,
          loss,
          penalty,
          gradient,
          working_set,
          x,
          x_centers,
          x_scales,
          y);
}

// Override for sparse matrices
void
Hybrid::run(Eigen::VectorXd& beta0,
            Eigen::MatrixXd& beta,
            Eigen::MatrixXd& eta,
            Clusters& clusters,
            const Eigen::ArrayXd& lambda,
            const std::unique_ptr<Loss>& loss,
            const SortedL1Norm& penalty,
            const Eigen::MatrixXd& gradient,
            const std::vector<int>& working_set,
            const Eigen::SparseMatrix<double>& x,
            const Eigen::VectorXd& x_centers,
            const Eigen::VectorXd& x_scales,
            const Eigen::MatrixXd& y)
{
  runImpl(beta0,
          beta,
          eta,
          clusters,
          lambda,
          loss,
          penalty,
          gradient,
          working_set,
          x,
          x_centers,
          x_scales,
          y);
}

} // namespace solvers
} // namespace slope
