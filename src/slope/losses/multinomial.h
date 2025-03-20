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
 * weights and working response. Assumes the response y is a one-hot encoded
 * matrix where each row sums to 1.
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
   * @param eta The predicted values (n x k matrix of linear predictors).
   * @param y The true labels (n x k matrix of one-hot encoded class
   * memberships).
   * @return The loss value.
   */
  double loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y);

  /**
   * @brief Calculates the dual for the multinomial loss function.
   * @param theta The dual variables (n x k matrix).
   * @param y The true labels (n x k matrix of one-hot encoded class
   * memberships).
   * @param w The weights vector.
   * @return The dual value.
   */
  double dual(const Eigen::MatrixXd& theta,
              const Eigen::MatrixXd& y,
              const Eigen::VectorXd& w);

  /**
   * @brief Calculates the residual for the multinomial loss function.
   * @param eta The predicted values (n x k matrix of linear predictors).
   * @param y The true labels (n x k matrix of one-hot encoded class
   * memberships).
   * @return The residual matrix (n x k).
   */
  Eigen::MatrixXd residual(const Eigen::MatrixXd& eta,
                           const Eigen::MatrixXd& y);

  /**
   * @brief Preprocesses the response for the Multinomial model
   * @param y Predicted values vector (n x 1) of integer class labels
   * @return Matrix of response (n x 1)
   */
  Eigen::MatrixXd preprocessResponse(const Eigen::MatrixXd& y);

  /**
   * @brief Updates the weights and working response for the multinomial
   * loss function. Currently not implemented since there is
   * no coordinate descent solver for the multinomial logistic regression
   * loss.
   * @param w The weights vector.
   * @param z The working response vector.
   * @param eta The predicted values (n x k matrix of linear predictors).
   * @param y The true labels (n x k matrix of one-hot encoded class
   * memberships).
   */
  void updateWeightsAndWorkingResponse(Eigen::VectorXd& w,
                                       Eigen::VectorXd& z,
                                       const Eigen::VectorXd& eta,
                                       const Eigen::VectorXd& y);

  /**
   * @brief The link function
   * @param mu Mean.
   * @return The result of applying the link function.
   */
  Eigen::MatrixXd link(const Eigen::MatrixXd& mu);

  /**
   * @brief The inverse link function, also known as the mean function.
   * @param eta
   * @return The softmax of the linear predictor: \f$
   * \frac{e^{\eta}}{\sum_{j=1}^{k} e^{\eta_j}} \f$
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
