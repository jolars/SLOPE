/**
 * @file
 * @brief Numerical solver class for SLOPE (Sorted L-One Penalized Estimation)
 */

#pragma once

#include "../jit_normalization.h"
#include "../losses/loss.h"
#include "../sorted_l1_norm.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>

namespace slope {

/**
 * @class SolverBase
 * @brief Abstract base class for SLOPE optimization solvers
 *
 * Provides the interface and common functionality for different SLOPE (Sorted
 * L-One Penalized Estimation) optimization algorithms. Derived classes
 * implement specific optimization strategies like coordinate descent or
 * proximal gradient descent.
 *
 * @see slope::PGD
 * @see slope::CD
 */
class SolverBase
{
public:
  /**
   * @brief Constructs a base solver for SLOPE optimization
   *
   * @param jit_normalization Type of just-in-time normalization to apply (None,
   * Center, Scale, or Both)
   * @param intercept Whether to fit an intercept term
   */
  SolverBase(JitNormalization jit_normalization, bool intercept)
    : jit_normalization(jit_normalization)
    , intercept(intercept)
  {
  }

  /// Default destructor
  virtual ~SolverBase() = default;

  /**
   * @brief Pure virtual function defining the solver's optimization routine
   *
   * @param beta0 Intercept terms for each response
   * @param beta Coefficients (size p x m)
   * @param eta Linear predictor matrix (n samples x m responses)
   * @param lambda Vector of regularization parameters
   * @param loss Pointer to loss function object
   * @param penalty Sorted L1 norm object for proximal operations
   * @param gradient Gradient matrix for loss function
   * @param working_set Vector of indices for active predictors
   * @param x Input feature matrix (n samples x p predictors)
   * @param x_centers Vector of feature means for centering
   * @param x_scales Vector of feature scales for normalization
   * @param y Response matrix (n samples x m responses)
   */
  virtual void run(Eigen::VectorXd& beta0,
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
                   const Eigen::MatrixXd& y) = 0;

  /**
   * @brief Pure virtual function defining the solver's optimization routine
   *
   * @param beta0 Intercept terms for each response
   * @param beta Coefficient vecttor (size p x m)
   * @param eta Linear predictor matrix (n samples x m responses)
   * @param lambda Vector of regularization parameters
   * @param loss Pointer to loss function object
   * @param penalty Sorted L1 norm object for proximal operations
   * @param gradient Gradient matrix for loss function
   * @param working_set Vector of indices for active predictors
   * @param x Input feature matrix (n samples x p predictors)
   * @param x_centers Vector of feature means for centering
   * @param x_scales Vector of feature scales for normalization
   * @param y Response matrix (n samples x m responses)
   */
  virtual void run(Eigen::VectorXd& beta0,
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
                   const Eigen::MatrixXd& y) = 0;

  /**
   * @brief Pure virtual function defining the solver's optimization routine
   *
   * @param beta0 Intercept terms for each response
   * @param beta Coefficient vecttor (size p x m)
   * @param eta Linear predictor matrix (n samples x m responses)
   * @param lambda Vector of regularization parameters
   * @param loss Pointer to loss function object
   * @param penalty Sorted L1 norm object for proximal operations
   * @param gradient Gradient matrix for loss function
   * @param working_set Vector of indices for active predictors
   * @param x Input feature matrix (n samples x p predictors)
   * @param x_centers Vector of feature means for centering
   * @param x_scales Vector of feature scales for normalization
   * @param y Response matrix (n samples x m responses)
   */
  virtual void run(Eigen::VectorXd& beta0,
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
                   const Eigen::MatrixXd& y) = 0;

  /**
   * @brief Pure virtual function defining the solver's optimization routine
   *
   * @param beta0 Intercept terms for each response
   * @param beta Coefficient vecttor (size p x m)
   * @param eta Linear predictor matrix (n samples x m responses)
   * @param lambda Vector of regularization parameters
   * @param loss Pointer to loss function object
   * @param penalty Sorted L1 norm object for proximal operations
   * @param gradient Gradient matrix for loss function
   * @param working_set Vector of indices for active predictors
   * @param x Input feature matrix (n samples x p predictors)
   * @param x_centers Vector of feature means for centering
   * @param x_scales Vector of feature scales for normalization
   * @param y Response matrix (n samples x m responses)
   */
  virtual void run(Eigen::VectorXd& beta0,
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
                   const Eigen::MatrixXd& y) = 0;

protected:
  JitNormalization jit_normalization; ///< JIT feature normalization strategy
  bool intercept;                     ///< If true, fits intercept term
};

} // namespace slope
