/**
 * @file
 * @brief Multinomial loss function implementation for SLOPE algorithm
 */

#pragma once

#include "loss.h"

namespace slope {
/**
 * @class Multinomial
 * @brief The Multinomial class represents a multinomial logistic regression
 * loss function.
 * @details The multinomial loss function is used for multi-class
 * classification problems. It calculates the loss, dual, residual, and updates
 * weights and working response. It uses the non-redundant formulation
 * of the loss with \f( K - 1\f) columns in the resulting
 * response matrix.
 */
class Multinomial : public Loss
{
public:
  explicit Multinomial()
    : Loss(1.0)
  {
  }
  /**
   * @brief Calculates the loss for the multinomial loss function.
   * @param eta The predicted values (n x m matrix of linear predictors).
   * @param y The true labels (n x m matrix of one-hot encoded class
   * memberships).
   * @return The loss value.
   */
  double loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y);

  /**
   * @brief Calculates the dual for the multinomial loss function.
   * @param theta The dual variables (n x m matrix).
   * @param y The true labels (n x m matrix).
   * @param w The weights vector.
   * @return The dual value.
   */
  double dual(const Eigen::MatrixXd& theta,
              const Eigen::MatrixXd& y,
              const Eigen::VectorXd& w);

  /**
   * @brief Preprocesses the response for the Multinomial model
   * @param y Vector of class labels (n x 1). Each entry is an integer
   * representing the class label from 0 to m.
   * @return Matrix of response (n x m)
   */
  Eigen::MatrixXd preprocessResponse(const Eigen::MatrixXd& y);

  /**
   * @brief Calculates the hessian diagonal
   *
   * @param eta Linear predictor
   * @param y Response
   * @return The hessian diagonal, a matrix of size (n x m)
   */
  Eigen::MatrixXd hessianDiagonal(const Eigen::MatrixXd& eta);

  /**
   * @brief The link function
   * @param mu Mean.
   * @return The result of applying the link function.
   */
  Eigen::MatrixXd link(const Eigen::MatrixXd& mu);

  /**
   * @brief The inverse link function
   * @param eta
   * @return The modified softmax of the linear predictor: \f$
   * \frac{e^{\eta}}{1 + \sum_{j=1}^{k} e^{\eta_j}} \f$
   */
  Eigen::MatrixXd inverseLink(const Eigen::MatrixXd& eta);

  /**
   * @brief Return predicted response, which is an integer class label based on
   *   the predicted probabilities.
   * @param eta The linear predictor
   * @return The predicted response
   */
  Eigen::MatrixXd predict(const Eigen::MatrixXd& eta);
};

} // namespace slope
