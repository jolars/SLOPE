/**
 * @file
 * @brief The declartion of the Objctive class and its subclasses, which
 * represent the data-fitting part of the composite loss function.
 */

#pragma once

#include <Eigen/Core>

namespace slope {

/**
 * Loss function interface
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
   * @brief Calculates the residual
   *
   * This function calculates the residual given the predicted values (eta) and
   * the true values (y).
   *
   * @param eta The predicted values.
   * @param y The true values.
   * @return The residual vector.
   */
  virtual Eigen::MatrixXd residual(const Eigen::MatrixXd& eta,
                                   const Eigen::MatrixXd& y) = 0;

  /**
   * @brief Updates the weights and working response
   *
   * This function updates the weights and working response given the predicted
   * values (eta) and the true values (y).
   *
   * @param w The weights to be updated.
   * @param z The working response to be updated.
   * @param eta The predicted values.
   * @param y The true values.
   */
  virtual void updateWeightsAndWorkingResponse(Eigen::VectorXd& w,
                                               Eigen::VectorXd& z,
                                               const Eigen::VectorXd& eta,
                                               const Eigen::VectorXd& y) = 0;

  /**
   * @brief Preprocess response
   *
   * @param y The response
   */
  virtual Eigen::MatrixXd preprocessResponse(const Eigen::MatrixXd& y) = 0;

  /**
   * @brief Updates the intercept with a
   * gradient descent update. Also updates the linear predictor (but not the
   * residual).
   * @param beta0 The current intercept
   * @param eta The current linear predictor
   * @param y The observed counts vector
   */
  virtual void updateIntercept(Eigen::VectorXd& beta0,
                               const Eigen::MatrixXd& eta,
                               const Eigen::MatrixXd& y)
  {
    Eigen::MatrixXd residual = this->residual(eta, y);
    beta0 -= residual.colwise().mean() / this->lipschitz_constant;
  };

  /**
   * @brief The link function
   * @param mu Mean.
   * @return The result of applying the link function.
   */
  virtual Eigen::MatrixXd link(const Eigen::MatrixXd& mu) = 0;

  /**
   * @brief Computes deviance, which is just twice the loss function
   * @param eta The predicted values.
   * @param y The true values.
   */
  virtual double deviance(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
  {
    return 2 * this->loss(eta, y);
  }

  /**
   * @brief Computes null deviance.
   * @param y The response matrix.
   * @param intercept Whether an intercept should be fit.
   */
  double nullDeviance(const Eigen::MatrixXd& y, const bool intercept)
  {
    int n = y.rows();
    int m = y.cols();

    Eigen::MatrixXd eta(n, m);

    if (intercept) {
      Eigen::RowVectorXd beta0(m);
      beta0 = this->link(y.colwise().mean());
      eta.rowwise() = beta0;
    } else {
      eta.setZero();
    }

    return deviance(eta, y);
  }

  /**
   * @brief Computes null deviance.
   * @param beta The response matrix.
   * @param intercept Whether an intercept should be fit.
   */
  virtual void adjustResponse(Eigen::MatrixXd&) {}

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
