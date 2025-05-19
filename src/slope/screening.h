/**
 * @file
 * @brief Screening rules for SLOPE regression optimization
 *
 * Implements feature screening methods to identify active and strong sets
 * of variables, reducing computational complexity in coordinate descent.
 */

#pragma once

#include "jit_normalization.h"
#include <Eigen/SparseCore>
#include <memory>
#include <vector>

namespace slope {

/**
 * @brief Identifies previously active variables
 *
 * @param beta Current coefficient matrix
 * @return std::vector<int> Indices of variables with non-zero coefficients
 */
std::vector<int>
activeSet(const Eigen::VectorXd& beta);

/**
 * @brief Determines the strong set using sequential strong rules
 *
 * @param gradient_prev Gradient from previous solution
 * @param lambda Current lambda sequence
 * @param lambda_prev Previous lambda sequence
 * @return std::vector<int> Indices of variables in the strong set
 */
std::vector<int>
strongSet(const Eigen::VectorXd& gradient_prev,
          const Eigen::ArrayXd& lambda,
          const Eigen::ArrayXd& lambda_prev);

/**
 * @class ScreeningRule
 * @brief Base class for screening rules in SLOPE.
 *
 * This abstract class defines the interface for screening rules to
 * accelerate optimization by reducing the problem size.
 */
class ScreeningRule
{
public:
  /**
   * @brief Virtual destructor
   */
  virtual ~ScreeningRule() = default;

  /**
   * @brief Initialize the screening rule at the start of the path algorithm.
   *
   * @param full_set The full set of feature indices
   * @param alpha_max_ind The index of the feature with maximum absolute
   * gradient
   * @return The initial working set
   */
  virtual std::vector<int> initialize(const std::vector<int>& full_set,
                                      int alpha_max_ind) = 0;

  /**
   * @brief Screen for the next path step.
   *
   * @param gradient The gradient vector
   * @param lambda_curr Current lambda values
   * @param lambda_prev Previous lambda values
   * @param beta Current beta coefficients
   * @param full_set Full set of features
   * @return Working set for the current path step
   */
  virtual std::vector<int> screen(Eigen::VectorXd& gradient,
                                  const Eigen::ArrayXd& lambda_curr,
                                  const Eigen::ArrayXd& lambda_prev,
                                  const Eigen::VectorXd& beta,
                                  const std::vector<int>& full_set) = 0;

  /**
   * @brief Check for KKT violations and update working set if necessary.
   *
   * @param gradient The gradient vector
   * @param beta Current beta coefficients
   * @param lambda_curr Current lambda values
   * @param working_set Current working set (will be updated if violations
   * found)
   * @param x Design matrix
   * @param residual Current residuals
   * @param x_centers Centers for normalization
   * @param x_scales Scales for normalization
   * @param jit_normalization Whether to use JIT normalization
   * @param full_set Full set of features
   * @return True if no violations found, false otherwise
   */
  virtual bool checkKktViolations(Eigen::VectorXd& gradient,
                                  const Eigen::VectorXd& beta,
                                  const Eigen::ArrayXd& lambda_curr,
                                  std::vector<int>& working_set,
                                  const Eigen::MatrixXd& x,
                                  const Eigen::MatrixXd& residual,
                                  const Eigen::VectorXd& x_centers,
                                  const Eigen::VectorXd& x_scales,
                                  JitNormalization jit_normalization,
                                  const std::vector<int>& full_set) = 0;
  /**
   * @brief Check for KKT violations with sparse matrix input
   * @param gradient The gradient vector
   * @param beta Current beta coefficients
   * @param lambda_curr Current lambda values
   * @param working_set Current working set (will be updated if violations
   * found)
   * @param x Design matrix (sparse format)
   * @param residual Current residuals
   * @param x_centers Centers for normalization
   * @param x_scales Scales for normalization
   * @param jit_normalization Whether to use JIT normalization
   * @param full_set Full set of features
   * @return True if no violations found, false otherwise
   */
  virtual bool checkKktViolations(Eigen::VectorXd& gradient,
                                  const Eigen::VectorXd& beta,
                                  const Eigen::ArrayXd& lambda_curr,
                                  std::vector<int>& working_set,
                                  const Eigen::SparseMatrix<double>& x,
                                  const Eigen::MatrixXd& residual,
                                  const Eigen::VectorXd& x_centers,
                                  const Eigen::VectorXd& x_scales,
                                  JitNormalization jit_normalization,
                                  const std::vector<int>& full_set) = 0;

  /**
   * @brief Check for KKT violations with sparse matrix input
   * @param gradient The gradient vector
   * @param beta Current beta coefficients
   * @param lambda_curr Current lambda values
   * @param working_set Current working set (will be updated if violations
   * found)
   * @param x Design matrix (sparse format)
   * @param residual Current residuals
   * @param x_centers Centers for normalization
   * @param x_scales Scales for normalization
   * @param jit_normalization Whether to use JIT normalization
   * @param full_set Full set of features
   * @return True if no violations found, false otherwise
   */
  virtual bool checkKktViolations(Eigen::VectorXd& gradient,
                                  const Eigen::VectorXd& beta,
                                  const Eigen::ArrayXd& lambda_curr,
                                  std::vector<int>& working_set,
                                  const Eigen::Map<Eigen::MatrixXd>& x,
                                  const Eigen::MatrixXd& residual,
                                  const Eigen::VectorXd& x_centers,
                                  const Eigen::VectorXd& x_scales,
                                  JitNormalization jit_normalization,
                                  const std::vector<int>& full_set) = 0;

  /**
   * @brief Check for KKT violations with sparse matrix input
   * @param gradient The gradient vector
   * @param beta Current beta coefficients
   * @param lambda_curr Current lambda values
   * @param working_set Current working set (will be updated if violations
   * found)
   * @param x Design matrix (sparse format)
   * @param residual Current residuals
   * @param x_centers Centers for normalization
   * @param x_scales Scales for normalization
   * @param jit_normalization Whether to use JIT normalization
   * @param full_set Full set of features
   * @return True if no violations found, false otherwise
   */
  virtual bool checkKktViolations(
    Eigen::VectorXd& gradient,
    const Eigen::VectorXd& beta,
    const Eigen::ArrayXd& lambda_curr,
    std::vector<int>& working_set,
    const Eigen::Map<Eigen::SparseMatrix<double>>& x,
    const Eigen::MatrixXd& residual,
    const Eigen::VectorXd& x_centers,
    const Eigen::VectorXd& x_scales,
    JitNormalization jit_normalization,
    const std::vector<int>& full_set) = 0;

  /**
   * @brief Get string representation of the screening rule
   * @return Name of the screening rule
   */
  virtual std::string toString() const = 0;

protected:
  /// Strong set of variables
  std::vector<int> strong_set;
};

/**
 * @class NoScreening
 * @brief No screening rule - uses all variables.
 */
class NoScreening : public ScreeningRule
{
public:
  std::vector<int> initialize(const std::vector<int>& full_set,
                              int alpha_max_ind) override;

  std::vector<int> screen(Eigen::VectorXd& gradient,
                          const Eigen::ArrayXd& lambda_curr,
                          const Eigen::ArrayXd& lambda_prev,
                          const Eigen::VectorXd& beta,
                          const std::vector<int>& full_set) override;

  bool checkKktViolations(Eigen::VectorXd& gradient,
                          const Eigen::VectorXd& beta,
                          const Eigen::ArrayXd& lambda_curr,
                          std::vector<int>& working_set,
                          const Eigen::MatrixXd& x,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          JitNormalization jit_normalization,
                          const std::vector<int>& full_set) override;

  bool checkKktViolations(Eigen::VectorXd& gradient,
                          const Eigen::VectorXd& beta,
                          const Eigen::ArrayXd& lambda_curr,
                          std::vector<int>& working_set,
                          const Eigen::SparseMatrix<double>& x,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          JitNormalization jit_normalization,
                          const std::vector<int>& full_set) override;

  bool checkKktViolations(Eigen::VectorXd& gradient,
                          const Eigen::VectorXd& beta,
                          const Eigen::ArrayXd& lambda_curr,
                          std::vector<int>& working_set,
                          const Eigen::Map<Eigen::MatrixXd>& x,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          JitNormalization jit_normalization,
                          const std::vector<int>& full_set) override;

  bool checkKktViolations(Eigen::VectorXd& gradient,
                          const Eigen::VectorXd& beta,
                          const Eigen::ArrayXd& lambda_curr,
                          std::vector<int>& working_set,
                          const Eigen::Map<Eigen::SparseMatrix<double>>& x,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          JitNormalization jit_normalization,
                          const std::vector<int>& full_set) override;

  std::string toString() const override;
};

/**
 * @class StrongScreening
 * @brief Implements strong screening rules for SLOPE.
 */
class StrongScreening : public ScreeningRule
{
public:
  std::vector<int> initialize(const std::vector<int>& full_set,
                              int alpha_max_ind) override;

  std::vector<int> screen(Eigen::VectorXd& gradient,
                          const Eigen::ArrayXd& lambda_curr,
                          const Eigen::ArrayXd& lambda_prev,
                          const Eigen::VectorXd& beta,
                          const std::vector<int>& full_set) override;

  bool checkKktViolations(Eigen::VectorXd& gradient,
                          const Eigen::VectorXd& beta,
                          const Eigen::ArrayXd& lambda_curr,
                          std::vector<int>& working_set,
                          const Eigen::MatrixXd& x,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          JitNormalization jit_normalization,
                          const std::vector<int>& full_set) override;

  bool checkKktViolations(Eigen::VectorXd& gradient,
                          const Eigen::VectorXd& beta,
                          const Eigen::ArrayXd& lambda_curr,
                          std::vector<int>& working_set,
                          const Eigen::Map<Eigen::MatrixXd>& x,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          JitNormalization jit_normalization,
                          const std::vector<int>& full_set) override;

  bool checkKktViolations(Eigen::VectorXd& gradient,
                          const Eigen::VectorXd& beta,
                          const Eigen::ArrayXd& lambda_curr,
                          std::vector<int>& working_set,
                          const Eigen::SparseMatrix<double>& x,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          JitNormalization jit_normalization,
                          const std::vector<int>& full_set) override;

  bool checkKktViolations(Eigen::VectorXd& gradient,
                          const Eigen::VectorXd& beta,
                          const Eigen::ArrayXd& lambda_curr,
                          std::vector<int>& working_set,
                          const Eigen::Map<Eigen::SparseMatrix<double>>& x,
                          const Eigen::MatrixXd& residual,
                          const Eigen::VectorXd& x_centers,
                          const Eigen::VectorXd& x_scales,
                          JitNormalization jit_normalization,
                          const std::vector<int>& full_set) override;

  std::string toString() const override;

private:
  template<typename MatrixType>
  bool checkKktViolationsImpl(Eigen::VectorXd& gradient,
                              const Eigen::VectorXd& beta,
                              const Eigen::ArrayXd& lambda_curr,
                              std::vector<int>& working_set,
                              const MatrixType& x,
                              const Eigen::MatrixXd& residual,
                              const Eigen::VectorXd& x_centers,
                              const Eigen::VectorXd& x_scales,
                              JitNormalization jit_normalization,
                              const std::vector<int>& full_set);
};

/**
 * @brief Creates a screening rule based on the provided type.
 *
 * @param screening_type Type of screening rule to create ("none" or "strong")
 * @return std::unique_ptr<ScreeningRule> A pointer to the created screening
 * rule
 */
std::unique_ptr<ScreeningRule>
createScreeningRule(const std::string& screening_type);

} // namespace slope
