/**
 * @file
 * @brief The declartion of the Objctive class and its subclasses, which
 * represent the data-fitting part of the composite loss function.
 */

#pragma once

#include <Eigen/Core>

namespace slope {

/**
 * @brief Loss function interface
 *
 * This class defines the interface for an loss function, which is used in
 * optimization algorithms. The loss function calculates the loss, dual,
 * residual, and updates the weights and working response.
 */
class Loss
{
public:
  /**
   * @brief Destructor for the Loss class.
   */
  virtual ~Loss() = default;

  /**
   * @brief Calculates the loss function
   *
   * This function calculates the loss function given the predicted values (eta)
   * and the true values (y).
   *
   * @param eta The predicted values.
   * @param y The true values.
   * @return The loss value.
   */
  virtual double loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y) = 0;

  /**
   * @brief Calculates the dual loss
   *
   * This function calculates the dual function given the estimated parameters
   * (theta) and the true values (y).
   *
   * @param theta The estimated parameters.
   * @param y The true values.
   * @param w Weights.
   * @return The dual value.
   */
  virtual double dual(const Eigen::MatrixXd& theta,
                      const Eigen::MatrixXd& y,
                      const Eigen::VectorXd& w) = 0;

  /**
   * @brief Calculates the generalized residual
   *
   * This function calculates the generalized residual given the linear
   * predictor (eta) and the true values (y). The generalized residual is the
   * same as the gradient of the loss function with respect to the linear
   * predictor (eta).
   *
   * @param eta Linear predictor.
   * @param y Response.
   * @return Residual.
   */
  Eigen::MatrixXd residual(const Eigen::MatrixXd& eta,
                           const Eigen::MatrixXd& y);
  // TODO: This function could be removed since it is just
  // g^{-1}(eta) - y for all families

  /**
   * @brief Calculates the hessian diagonal of the loss function
   *
   * This function calculates the hessian diagonal, with respect to
   * the linear predictor eta. \f(\frac{\partial^2 L}{\partial \eta^2}\f)
   *
   * @param eta Linear predictor.
   * @param y Response.
   * @return The residual.
   */
  virtual Eigen::MatrixXd hessianDiagonal(const Eigen::MatrixXd& eta) = 0;

  /**
   * @brief Updates weights and working response
   *
   * This function updates the weights and working response of the
   * quadratic expansion of the loss function.
   *
   * @param w Working weights.
   * @param z Working response.
   * @param eta Linear predictor.
   * @param y Response.
   */
  virtual void updateWeightsAndWorkingResponse(Eigen::MatrixXd& w,
                                               Eigen::MatrixXd& z,
                                               const Eigen::MatrixXd& eta,
                                               const Eigen::MatrixXd& y);

  /**
   * @brief Preprocess response
   *
   * @param y The response
   */
  virtual Eigen::MatrixXd preprocessResponse(const Eigen::MatrixXd& y) = 0;

  /**
   * @brief Updates the intercept with a
   * gradient descent update.
   *
   * @param beta0 The current intercept.
   * @param eta The current linear predictor.
   * @param y The observed counts vector.
   */
  virtual void updateIntercept(Eigen::VectorXd& beta0,
                               const Eigen::MatrixXd& eta,
                               const Eigen::MatrixXd& y)
  {
    Eigen::MatrixXd residual = this->residual(eta, y);
    beta0 -= residual.colwise().mean() / this->lipschitz_constant;
  };
  // TODO: This function is not currently used anywhere and could
  // probably be removed.

  /**
   * @brief The link function.
   * @param mu Mean.
   * @return The result of applying the link function, which depends on
   *   the loss function used.
   */
  virtual Eigen::MatrixXd link(const Eigen::MatrixXd& mu) = 0;

  /**
   * @brief The inverse link function, also known as the mean function.
   * @param eta Mean
   * @return The result of applying the inverse link function, which depends on
   *   the loss function used.
   */
  virtual Eigen::MatrixXd inverseLink(const Eigen::MatrixXd& eta) = 0;

  /**
   * @brief Return predicted response
   * @param eta The linear predictor.
   * @return The predicted response.
   */
  virtual Eigen::MatrixXd predict(const Eigen::MatrixXd& eta) = 0;

  /**
   * @brief Computes deviance, which is 2 times the difference between the
   * loglikelihood of the model and the loglikelihood of the null
   * (intercept-only) model.
   * @param eta The predicted values.
   * @param y The true values.
   */
  virtual double deviance(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
  {
    return 2.0 * (loss(eta, y) - loss(link(y), y));
  }

protected:
  /**
   * @brief Constructs an loss function with specified Lipschitz constant
   * @param lipschitz_constant The Lipschitz constant for the loss function
   *
   * The Lipschitz constant is used to ensure convergence in gradient-based
   * optimization by bounding the rate of change of the gradient.
   */
  explicit Loss(double lipschitz_constant)
    : lipschitz_constant(lipschitz_constant)
  {
  }

private:
  const double lipschitz_constant;
};

} // namespace slope
