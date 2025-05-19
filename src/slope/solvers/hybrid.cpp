/**
 * @file
 * @brief Hybrid solver implementation for SLOPE
 */

#include "hybrid.h"
#include "../losses/loss.h"
#include "../sorted_l1_norm.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <memory>

namespace slope {

// Override for dense matrices
void
Hybrid::run(Eigen::VectorXd& beta0,
            Eigen::VectorXd& beta,
            Eigen::MatrixXd& eta,
            const Eigen::ArrayXd& lambda,
            const std::unique_ptr<Loss>& loss,
            const SortedL1Norm& penalty,
            const Eigen::VectorXd& gradient,
            const std::vector<int>& working_set,
            const Eigen::MatrixXd& x,
            const Eigen::VectorXd& x_centers,
            const Eigen::VectorXd& x_scales,
            const Eigen::MatrixXd& y)
{
  runImpl(beta0,
          beta,
          eta,
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
            Eigen::VectorXd& beta,
            Eigen::MatrixXd& eta,
            const Eigen::ArrayXd& lambda,
            const std::unique_ptr<Loss>& loss,
            const SortedL1Norm& penalty,
            const Eigen::VectorXd& gradient,
            const std::vector<int>& working_set,
            const Eigen::SparseMatrix<double>& x,
            const Eigen::VectorXd& x_centers,
            const Eigen::VectorXd& x_scales,
            const Eigen::MatrixXd& y)
{
  runImpl(beta0,
          beta,
          eta,
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

void
Hybrid::run(Eigen::VectorXd& beta0,
            Eigen::VectorXd& beta,
            Eigen::MatrixXd& eta,
            const Eigen::ArrayXd& lambda,
            const std::unique_ptr<Loss>& loss,
            const SortedL1Norm& penalty,
            const Eigen::VectorXd& gradient,
            const std::vector<int>& working_set,
            const Eigen::Map<Eigen::MatrixXd>& x,
            const Eigen::VectorXd& x_centers,
            const Eigen::VectorXd& x_scales,
            const Eigen::MatrixXd& y)
{
  runImpl(beta0,
          beta,
          eta,
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

void
Hybrid::run(Eigen::VectorXd& beta0,
            Eigen::VectorXd& beta,
            Eigen::MatrixXd& eta,
            const Eigen::ArrayXd& lambda,
            const std::unique_ptr<Loss>& loss,
            const SortedL1Norm& penalty,
            const Eigen::VectorXd& gradient,
            const std::vector<int>& working_set,
            const Eigen::Map<Eigen::SparseMatrix<double>>& x,
            const Eigen::VectorXd& x_centers,
            const Eigen::VectorXd& x_scales,
            const Eigen::MatrixXd& y)
{
  runImpl(beta0,
          beta,
          eta,
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

} // namespace slope
