/**
 * @file
 * @brief Quadratic loss function implementation for SLOPE algorithm
 * @details This file contains the Quadratic class which implements a
 * Quadratic loss function used in the SLOPE (Sorted L-One Penalized
 * Estimation) algorithm. The Quadratic loss function is particularly useful for
 * regression problems with normally distributed errors.
 */

#pragma once

#include "loss.h"

namespace slope {

/**
 * @class Quadratic
 * @brief Implementation of the Quadratic loss function
 * @details The Quadratic class provides methods for computing loss, dual
 * function, residuals, and weight updates for the Quadratic case in the SLOPE
 * algorithm. It is particularly suited for regression problems where the error
 * terms are assumed to follow a normal distribution.
 *
 * @note This class inherits from the base Loss class and implements
 * all required virtual functions.
 */
class Quadratic : public Loss
{
public:
  explicit Quadratic()
    : Loss(1.00)
  {
  }
  /**
   * @brief Calculates the quadratic (least-squares) loss.
   * @details Computes the squared error loss between predicted and actual
   * values, normalized by twice the number of observations.
   *
   * @param eta Vector of predicted values (n x 1)
   * @param y Matrix of actual values (n x 1)
   * @return Double precision loss value
   *
   * @note The loss is calculated as: \f$ \frac{1}{2n} \sum_{i=1}^n (\eta_i -
   * y_i)^2 \f$
   */
  double loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y);

  /**
   * @brief Computes the dual function for the quadratic loss
   * @details Calculates the Fenchel conjugate of the quadratic loss function
   *
   * @param theta Dual variables vector (n x 1)
   * @param y Observed values vector (n x 1)
   * @param w Observation weights vector (n x 1)
   * @return Double precision dual value
   *
   * @see loss() for the primal function
   */
  double dual(const Eigen::MatrixXd& theta,
              const Eigen::MatrixXd& y,
              const Eigen::VectorXd& w);

  /**
   * @brief Calculates hessian diagonal
   *
   * @param eta Linear predictor
   * @param y Response
   * @return A matrix of ones (n x m)
   */
  Eigen::MatrixXd hessianDiagonal(const Eigen::MatrixXd& eta);

  /**
   * @brief Preprocesses the response for the quadratic model
   * @details Doesn't perform any transformation on the response.
   *
   * @param y Responnse vector (n x 1)
   * @return Modified response
   */
  Eigen::MatrixXd preprocessResponse(const Eigen::MatrixXd& y);

  /**
   * @brief Updates weights and working response for IRLS algorithm
   * @details For quadratic case, weights are set to 1 and working response
   * equals the original response. This implementation is particularly simple
   * compared to other GLM families.
   *
   * @param[out] w Weights vector to be updated (n x 1)
   * @param[out] z Working response vector to be updated (n x 1)
   * @param[in] eta Current predictions vector (n x 1)
   * @param[in] y Matrix of observed values (n x 1)
   *
   * @note For quadratic regression, this is particularly simple as weights
   * remain constant and working response equals the original response
   */
  void updateWeightsAndWorkingResponse(Eigen::MatrixXd& w,
                                       Eigen::MatrixXd& z,
                                       const Eigen::MatrixXd& eta,
                                       const Eigen::MatrixXd& y);

  /**
   * @brief The link function
   * @param mu Mean of the distribution.
   * @return The identity function.
   */
  Eigen::MatrixXd link(const Eigen::MatrixXd& mu);

  /**
   * @brief The link function, also known as the mean function.
   * @param eta Linear predictor.
   * @return The identity function.
   */
  Eigen::MatrixXd inverseLink(const Eigen::MatrixXd& eta);

  /**
   * @brief Return predicted response, which is the same as the linear predictor
   * @param eta The linear predictor
   * @return The predicted response.
   */
  Eigen::MatrixXd predict(const Eigen::MatrixXd& eta);
};

} // namespace slope
