/**
 * @file
 * @brief Hybrid numerical solver for SLOPE combining coordinate descent and
 * proximal gradient descent
 */

#pragma once

#include "../clusters.h"
#include "../losses/loss.h"
#include "../sorted_l1_norm.h"
#include "hybrid_cd.h"
#include "pgd.h"
#include "solver.h"
#include <memory>

namespace slope {

/**
 * @brief Hybrid CD-PGD solver for SLOPE
 *
 * This solver alternates between coordinate descent (CD) and proximal gradient
 * descent (PGD) steps to solve the SLOPE optimization problem. Because
 * the SLOPE problem is non-separable, CD steps do not work out-of-the-box.
 * In particular, they cannot be used to split he clusters (or would need
 * to be augmented with additional steps). Instead, we use a hybrid method:
 * - PGD: Splits (or merges) clusters
 * - CD: Coordinate descent over the clusters, with good performance but
 *   can only merge, not split, clusters.
 * cases
 *
 * The switching between methods is controlled by the cd_iterations parameter,
 * which determines how often PGD steps are taken versus CD steps.
 */
class Hybrid : public SolverBase
{
public:
  /**
   * @brief Constructs Hybrid solver for SLOPE optimization
   * @param jit_normalization Feature normalization strategy
   * @param intercept If true, fits intercept term
   * @param update_clusters If true, updates clusters during optimization
   * @param cd_iterations Frequency of proximal gradient descent updates
   */
  Hybrid(JitNormalization jit_normalization,
         bool intercept,
         bool update_clusters,
         int cd_iterations)
    : SolverBase(jit_normalization, intercept)
    , update_clusters(update_clusters)
    , cd_iterations(cd_iterations)
  {
  }

  /// @copydoc SolverBase::run
  void run(Eigen::VectorXd& beta0,
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
           const Eigen::MatrixXd& y) override;

  /// @copydoc SolverBase::run
  void run(Eigen::VectorXd& beta0,
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
           const Eigen::MatrixXd& y) override;

  /// @copydoc SolverBase::run
  void run(Eigen::VectorXd& beta0,
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
           const Eigen::MatrixXd& y) override;

  /// @copydoc SolverBase::run
  void run(Eigen::VectorXd& beta0,
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
           const Eigen::MatrixXd& y) override;

private:
  /**
   * @brief Implementation of the hybrid solver algorithm
   *
   * @tparam MatrixType Type of the design matrix
   * @param beta0 Intercept term (scalar)
   * @param beta Coefficients
   * @param eta Linear predictor
   * @param loss Pointer to the loss function
   * @param penalty SLOPE penalty object
   * @param x Design matrix
   * @param x_centers Feature centers for standardization
   * @param x_scales Feature scales for standardization
   * @param y Response variable
   */
  template<typename MatrixType>
  void runImpl(Eigen::VectorXd& beta0,
               Eigen::VectorXd& beta,
               Eigen::MatrixXd& eta,
               const Eigen::ArrayXd& lambda,
               const std::unique_ptr<Loss>& loss,
               const SortedL1Norm& penalty,
               const Eigen::VectorXd& gradient_in,
               const std::vector<int>& working_set,
               const MatrixType& x,
               const Eigen::VectorXd& x_centers,
               const Eigen::VectorXd& x_scales,
               const Eigen::MatrixXd& y)
  {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int n = x.rows();
    const int m = eta.cols();

    PGD pgd_solver(jit_normalization, intercept, "pgd");

    // Run proximal gradient descent
    pgd_solver.run(beta0,
                   beta,
                   eta,
                   lambda,
                   loss,
                   penalty,
                   gradient_in,
                   working_set,
                   x,
                   x_centers,
                   x_scales,
                   y);

    Clusters clusters(beta);

    // TODO: Make these parameters and initialize once
    MatrixXd w = MatrixXd::Ones(n, m);
    MatrixXd z = y;

    loss->updateWeightsAndWorkingResponse(w, z, eta, y);

    MatrixXd residual = eta - z;

    for (int it = 0; it < this->cd_iterations; ++it) {
      double old_obj =
        computeObjective(penalty, beta, residual, w, lambda, working_set);

      // Store old values to revert if no progress is made
      Clusters old_clusters = clusters;
      Eigen::MatrixXd old_residual = residual;
      Eigen::VectorXd old_beta = beta;
      Eigen::VectorXd old_beta0 = beta0;

      coordinateDescent(beta0,
                        beta,
                        residual,
                        clusters,
                        lambda,
                        x,
                        w,
                        x_centers,
                        x_scales,
                        this->intercept,
                        this->jit_normalization,
                        this->update_clusters);

      double new_obj =
        computeObjective(penalty, beta, residual, w, lambda, working_set);

      if (!std::isfinite(new_obj) || new_obj > old_obj) {
        // No progress, revert to previous state
        clusters = old_clusters;
        residual = old_residual;
        beta = old_beta;
        beta0 = old_beta0;

        break;
      }
    }

    // The residual is kept up to date, but not eta. So we need to compute
    // it here.
    eta = residual + z;
    // TODO: register convergence status
  }

  double computeObjective(const SortedL1Norm& penalty,
                          const Eigen::VectorXd& beta,
                          const Eigen::MatrixXd& residual,
                          const Eigen::MatrixXd& w,
                          const Eigen::ArrayXd& lambda,
                          const std::vector<int>& working_set)
  {
    double val =
      0.5 * (residual.array().square() * w.array()).sum() / residual.rows() +
      penalty.eval(beta(working_set), lambda.head(working_set.size()));

    return val;
  }

  // TODO: These should be used in the PGD solver and taken as arguments to the
  // Hybrid solver and not just set and ignored here.
  double pgd_learning_rate =
    1.0; ///< Learning rate for proximal gradient descent steps
  double pgd_learning_rate_decr =
    0.5; ///< Learning rate decrease factor on failed PGD steps

  bool update_clusters = false; ///< If true, updates clusters during CD steps
  int cd_iterations = 10;       ///< Number of CD iterations per hybrid step
};

} // namespace slope
