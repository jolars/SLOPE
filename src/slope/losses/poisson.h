/**
 * @file
 * @brief Poisson loss function implementation for SLOPE algorithm
 * @details This file contains the Poisson class which implements a
 * Poisson loss function used in the SLOPE (Sorted L-One Penalized Estimation)
 * algorithm.
 */

#pragma once

#include "loss.h"

namespace slope {
/**
 * @class Poisson
 * @brief The Poisson class represents a Poisson regression loss function.
 * @details The Poisson regression loss function is used for modeling count
 * data. It assumes the response variable follows a Poisson distribution. The
 * log-likelihood for a single observation is:
 * \f[ \ell(y_i|\eta_i) = y_i\eta_i - e^{\eta_i} - \log(y_i!) \f]
 * where \f$\eta_i\f$ is the linear predictor and \f$y_i\f$ is the observed
 * count.
 */
class Poisson : public Loss
{

public:
  explicit Poisson()
    : Loss(std::numeric_limits<double>::infinity())
  {
  }

  /**
   * @brief Calculates the negative log-likelihood loss for the Poisson
   * regression.
   * @param eta The linear predictor vector \f$\eta\f$
   * @param y The observed counts vector \f$y\f$
   * @return The negative log-likelihood: \f$-\sum_i(y_i\eta_i - e^{\eta_i})\f$
   * (ignoring constant terms)
   */
  double loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y) override;

  /**
   * @brief Calculates the Fenchel conjugate (dual) of the Poisson loss.
   * @param theta The dual variables
   * @param y The observed counts vector
   * @param w The weights vector
   * @return The dual objective value
   */
  double dual(const Eigen::MatrixXd& theta,
              const Eigen::MatrixXd& y,
              const Eigen::VectorXd& w) override;

  /**
   * @brief Calculates hessian diagonal
   *
   * @param eta Linear predictor
   * @param y Response
   * @return A matrix of ones (n x m)
   */
  Eigen::MatrixXd hessianDiagonal(const Eigen::MatrixXd& eta) override;

  /**
   * @brief Preprocesses the response for the Poisson model
   * @details Checks if the response is non-negative and throws an error
   * otherwise.
   *
   * @param y Predicted values vector (n x 1)
   * @return Vector of residuals (n x 1)
   */
  Eigen::MatrixXd preprocessResponse(const Eigen::MatrixXd& y) override;

  /**
   * @brief Updates the intercept with a
   * gradient descent update. Unlike the Quadratic and Logistic cases,
   * the Poisson regression intercept update is not quite as simple since the
   * gradient is not Lipschitz continuous. Instead we use
   * a backtracking line search here.
   * @param beta0 The current intercept
   * @param eta The current linear predictor
   * @param y The observed counts vector
   */
  void updateIntercept(Eigen::VectorXd& beta0,
                       const Eigen::MatrixXd& eta,
                       const Eigen::MatrixXd& y) override;

  /**
   * @brief The link function
   * @param mu Mean.
   * @return \f$ \log(\mu) \f$
   */
  Eigen::MatrixXd link(const Eigen::MatrixXd& mu) override;

  /**
   * @brief The inverse link function, also known as the mean function.
   * @param eta Linear predictor
   * @return \f$ \exp(\eta) \f$
   */
  Eigen::MatrixXd inverseLink(const Eigen::MatrixXd& eta) override;

  /**
   * @brief Return predicted response, that is 0 or 1 depending on
   *   the predicted probabilities.
   * @param eta The linear predictor
   * @return The predicted response, which is the same as the inverse link
   *   in the case of the Poisson loss.
   */
  Eigen::MatrixXd predict(const Eigen::MatrixXd& eta) override;
};

} // namespace slope
