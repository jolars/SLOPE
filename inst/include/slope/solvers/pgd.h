/**
 * @file
 * @brief Proximal Gradient Descent solver implementation for SLOPE
 */

#pragma once

#include "../losses/loss.h"
#include "../math.h"
#include "../sorted_l1_norm.h"
#include "solver.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <memory>

namespace slope {

/**
 * @brief Proximal Gradient Descent solver for SLOPE optimization
 *
 * This solver implements the proximal gradient descent algorithm with line
 * search for solving the SLOPE optimization problem. It uses backtracking line
 * search to automatically adjust the learning rate for optimal convergence.
 */
class PGD : public SolverBase
{
public:
  /**
   * @brief Constructs Proximal Gradient Descent solver for SLOPE optimization
   * @param jit_normalization Feature normalization strategy
   * @param intercept If true, fits intercept term
   * @param update_type Type of update strategy to use
   */
  PGD(JitNormalization jit_normalization,
      bool intercept,
      const std::string& update_type)
    : SolverBase(jit_normalization, intercept)
    , learning_rate(1.0)
    , learning_rate_decr(0.5)
    , update_type{ update_type }
    , t(1.0)
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
  template<typename MatrixType>
  void runImpl(Eigen::VectorXd& beta0,
               Eigen::VectorXd& beta,
               Eigen::MatrixXd& eta,
               const Eigen::ArrayXd& lambda,
               const std::unique_ptr<Loss>& loss,
               const SortedL1Norm& penalty,
               const Eigen::VectorXd& gradient,
               const std::vector<int>& working_set,
               const MatrixType& x,
               const Eigen::VectorXd& x_centers,
               const Eigen::VectorXd& x_scales,
               const Eigen::MatrixXd& y)
  {
    using Eigen::all;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    int n_working = working_set.size();

    Eigen::VectorXd beta_old = beta(working_set);
    Eigen::VectorXd beta_diff(n_working);
    Eigen::VectorXd beta0_old = beta0;

    double g_old = loss->loss(eta, y);
    double t_old = t;

    Eigen::MatrixXd residual = loss->residual(eta, y);
    Eigen::VectorXd intercept_grad = residual.colwise().mean();
    Eigen::VectorXd grad_working = gradient(working_set);

    int line_search_iter = 0;
    const int max_line_search_iter = 100; // maximum iterations before exit

    while (true) {
      if (intercept) {
        beta0 = beta0_old - this->learning_rate * intercept_grad;
      }

      beta(working_set) =
        penalty.prox(beta_old - this->learning_rate * grad_working,
                     this->learning_rate * lambda.head(n_working));

      beta_diff = beta(working_set) - beta_old;
      double beta_diff_norm =
        beta_diff.dot(grad_working) +
        (1.0 / (2 * this->learning_rate)) * beta_diff.squaredNorm();

      if (intercept) {
        Eigen::VectorXd beta0_diff = beta0 - beta0_old;

        beta_diff_norm +=
          intercept_grad.dot(beta0_diff) +
          (1.0 / (2 * this->learning_rate)) * beta0_diff.squaredNorm();
      }

      eta = linearPredictor(x,
                            working_set,
                            beta0,
                            beta,
                            x_centers,
                            x_scales,
                            jit_normalization,
                            intercept);

      double g = loss->loss(eta, y);
      double q = g_old + beta_diff_norm;

      if (q >= g * (1 - 1e-12) || this->learning_rate < 1e-12) {
        this->learning_rate *= 1.1;
        break;
      } else {
        this->learning_rate *= this->learning_rate_decr;
      }

      line_search_iter++;
      if (line_search_iter >= max_line_search_iter) {
        break;
      }
    }

    if (update_type == "fista") {
      if (beta_prev.size() == 0) {
        beta_prev = beta;
      }

      this->t = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t_old * t_old));
      Eigen::VectorXd beta_current = beta(working_set);

      beta(working_set) +=
        (beta_current - beta_prev(working_set)) * (t_old - 1.0) / this->t;
      beta_prev(working_set) = beta_current;

      eta = linearPredictor(x,
                            working_set,
                            beta0,
                            beta,
                            x_centers,
                            x_scales,
                            jit_normalization,
                            intercept);
    }
  }

  double learning_rate;      ///< Current learning rate for gradient steps
  double learning_rate_decr; ///< Learning rate decrease factor for line search
  std::string update_type;   ///< Update type for PGD
  double t;                  ///< FISTA step size
  Eigen::VectorXd beta_prev; ///< Old beta values
};

} // namespace slope
