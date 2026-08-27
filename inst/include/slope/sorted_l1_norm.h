/**
 * @file
 * @brief The declaration of the SortedL1Norm class
 */

#pragma once

#include <Eigen/Core>

namespace slope {

/**
 * @brief Class representing the Sorted L1 Norm.
 */
class SortedL1Norm
{
public:
  /**
   * @brief Evaluates the Sorted L1 Norm.
   * @param beta The beta parameter.
   * @param lambda The regularization weights.
   * @return The evaluation result.
   */
  double eval(const Eigen::VectorXd& beta, const Eigen::ArrayXd& lambda) const;

  /**
   * @brief Computes the proximal operator of the Sorted L1 Norm.
   * @param beta The beta parameter.
   * @param lambda The regulariation weights.
   * @return The proximal operator result.
   */
  Eigen::MatrixXd prox(const Eigen::VectorXd& beta,
                       const Eigen::ArrayXd& lambda) const;

  /**
   * @brief Computes the dual norm of a vector.
   * @param a The vector.
   * @param lambda The regularization weights.
   * @return The dual norm.
   */
  double dualNorm(const Eigen::VectorXd& a, const Eigen::ArrayXd& lambda) const;

  /**
   * @brief Computes the dual norm with a numerical denominator floor.
   *
   * This overload is used for restricted working-set convergence checks, where
   * treating a tiny positive penalty as exactly zero can prevent the working
   * set from expanding. Full dual values use the unfloored overload.
   *
   * @param a The vector whose dual norm is computed.
   * @param lambda The regularization weights.
   * @param denominator_floor Minimum cumulative penalty denominator.
   * @return The numerically regularized dual norm.
   */
  double dualNorm(const Eigen::VectorXd& a,
                  const Eigen::ArrayXd& lambda,
                  double denominator_floor) const;
};

} // namespace slope
