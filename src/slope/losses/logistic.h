/**
 * @file
 * @brief Logistic loss function implementation for SLOPE algorithm
 */

#pragma once

#include "loss.h"

namespace slope {
/**
 * @class Logistic
 * @brief The Logistic class represents a logistic loss function.
 * @details The logistic loss function is used for binary classification
 * problems. It calculates the loss, dual, residual, and updates weights and
 * working response.
 */
class Logistic : public Loss
{
public:
  explicit Logistic()
    : Loss(0.25)
  {
  }

  /**
   * @brief Calculates the loss for the logistic loss function.
   * @param eta The predicted values.
   * @param y The true labels.
   * @return The loss value.
   */
  double loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y);

  /**
   * @brief Calculates the dual for the logistic loss function.
   * @param theta The dual variables.
   * @param y The true labels.
   * @param w Weights
   * @return The dual value.
   */
  double dual(const Eigen::MatrixXd& theta,
              const Eigen::MatrixXd& y,
              const Eigen::VectorXd& w);

  /**
   * @brief Calculates the residual for the logistic loss function.
   * @param eta The predicted values.
   * @param y The true labels.
   * @return The residual vector.
   */
  Eigen::MatrixXd residual(const Eigen::MatrixXd& eta,
                           const Eigen::MatrixXd& y);

  /**
   * @brief Preprocesses the response for the quadratic model
   * @details Checks if the response is in {0, 1} and converts it otherwise
   *
   * @param y Response vector (in {0,1})
   * @return Modified response.
   */
  Eigen::MatrixXd preprocessResponse(const Eigen::MatrixXd& y);

  /**
   * @brief Updates the weights and working response for the logistic loss
   * function.
   * @param w The weights.
   * @param z The working response.
   * @param eta The predicted values.
   * @param y The true labels.
   */
  void updateWeightsAndWorkingResponse(Eigen::VectorXd& w,
                                       Eigen::VectorXd& z,
                                       const Eigen::VectorXd& eta,
                                       const Eigen::VectorXd& y);

  /**
   * @brief The link function
   * @param mu Mean
   * @return \f$ \log(\frac{\mu}{1 - \mu}) \f$
   */
  Eigen::MatrixXd link(const Eigen::MatrixXd& mu);

  /**
   * @brief The inverse link function, also known as the mean function.
   * @param eta Mean
   * @return \f$ \frac{1}{1 + \exp(-\eta)} \f$
   */
  Eigen::MatrixXd inverseLink(const Eigen::MatrixXd& eta);

  /**
   * @brief Return predicted response, that is 0 or 1 depending on
   *   the predicted probabilities.
   * @param eta The linear predictor
   * @return The predicted response
   */
  Eigen::MatrixXd predict(const Eigen::MatrixXd& eta);
};

} // namespace slope
