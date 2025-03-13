/*
 * @file
 * @brief The actual function that fits SLOPE
 */

#pragma once

#include "slope_fit.h"
#include "slope_path.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cassert>
#include <optional>

/** @namespace slope
 *  @brief Namespace containing SLOPE regression implementation
 */
namespace slope {

/**
 * Class representing SLOPE (Sorted L-One Penalized Estimation) optimization.
 *
 * This class implements the SLOPE algorithm for regularized regression
 * problems. It supports different loss functions (quadratic, logistic, poisson)
 * and provides functionality for fitting models with sorted L1 regularization
 * along a path of regularization parameters.
 */
class Slope
{
public:
  /**
   * Default constructor for the Slope class.
   *
   * Initializes the Slope object with default parameter values.
   */
  Slope()
    : intercept(true)
    , modify_x(false)
    , update_clusters(false)
    , collect_diagnostics(false)
    , return_clusters(false)
    , alpha_min_ratio(-1) // TODO: Use std::optional for alpha_min_ratio
    , dev_change_tol(1e-5)
    , dev_ratio_tol(0.999)
    , learning_rate_decr(0.5)
    , q(0.1)
    , tol(1e-4)
    , max_it(1e4)
    , path_length(100)
    , cd_iterations(10)
    , max_clusters(std::optional<int>())
    , lambda_type("bh")
    , centering_type("mean")
    , scaling_type("sd")
    , loss_type("quadratic")
    , screening_type("strong")
    , solver_type("auto")
  {
  }

  /**
   * @brief Sets the numerical solver used to fit the model.
   *
   * @param solver One of "auto", "pgd", "fista", or "hybrid". In the first case
   * (the default), the solver is automatically selected based on availability
   * of the hybrid solver, which currently means that the hybrid solver is used
   * everywhere except for the multinomial loss.
   */
  void setSolver(const std::string& solver);

  /**
   * @brief Sets the intercept flag.
   *
   * @param intercept Should an intercept be fitted?
   */
  void setIntercept(bool intercept);

  /**
   * @brief Sets normalization type for the design matrix.
   *
   * @param type Type of normalization: one of "standardization" or "none".
   * @see setCentering() For setting centers, specifically
   * @see setScaling() For setting scales, specifically
   */
  void setNormalization(const std::string& type);

  /**
   * @brief Sets the update clusters flag.
   *
   * @param update_clusters Selects whether the coordinate descent keeps the
   * clusters updated.
   */
  void setUpdateClusters(bool update_clusters);

  /**
   * @brief Sets the update clusters flag.
   *
   * @param return_clusters Selects whether the fitted model should return
   * cluster information.
   */
  void setReturnClusters(const bool return_clusters);

  /**
   * @brief Sets the alpha min ratio.
   *
   * @param alpha_min_ratio The value to set for the alpha min ratio. A negative
   * value means that the program automatically chooses 1e-4 if the number of
   * observations is larger than the number of features and 1e-2 otherwise.
   */
  void setAlphaMinRatio(double alpha_min_ratio);

  /**
   * @brief Sets the learning rate decrement.
   *
   * @param learning_rate_decr The value to set for the learning rate decrement
   * for the proximal gradient descent step.
   */
  void setLearningRateDecr(double learning_rate_decr);

  /**
   * @brief Sets the q value.
   *
   * @param q The value to set for the q value for use in automatically
   * generating the lambda sequence. values between 0 and 1 are allowed..
   */
  void setQ(double q);

  /**
   * @brief Sets OSCAR parameters.
   *
   * @param theta1 Parameter for OSCAR.
   * @param theta2 Parameter for OSCAR.
   */
  void setOscarParameters(const double theta1, const double theta2);

  /**
   * @brief Sets the tolerance value.
   *
   * @param tol The value to set for the tolerance value. Must be positive.
   */
  void setTol(double tol);

  /**
   * @brief Sets the maximum number of iterations.
   *
   * @param max_it The value to set for the maximum number of iterations. Must
   * be positive. If negative (the default), then the value will be decided by
   * the solver.
   */
  void setMaxIterations(int max_it);

  /**
   * @brief Sets the path length.
   *
   * @param path_length The value to set for the path length.
   */
  void setPathLength(int path_length);

  /**
   * @brief Sets the frequence of proximal gradient descent steps.
   *
   * @param cd_iterations Number of inner coordinate descent iterations to
   * perform in each iteration of the hybrid solver. If set to 0, the solver
   * will resolve into pure PGD.
   */
  void setHybridCdIterations(int cd_iterations);

  /**
   * @brief Sets the lambda type for regularization weights.
   *
   * @param lambda_type The method used to compute regularization weights.
   * Currently "bh" (Benjamini-Hochberg), "quadratic", "oscar", and "lasso" are
   * supported.
   */
  void setLambdaType(const std::string& lambda_type);

  /**
   * @brief Sets the loss function type.
   *
   * @param loss_type The type of loss function to use. Supported values
   * are:
   *                 - "quadratic": Quadratic regression
   *                 - "logistic": Logistic regression
   *                 - "poisson": Poisson regression
   *                 - "multinomial": Multinomial logistic regression
   */
  void setLoss(const std::string& loss_type);

  /**
   * @brief Sets the type of feature screening used, which discards predictors
   * that are unlikely to be active.
   *
   * @param screening_type Type of screening. Supported values are:
   * are:
   *   - "strong": Strong screening rule ()
   *   - "none": No screening
   */
  void setScreening(const std::string& screening_type);

  /**
   * @brief Controls if `x` should be modified-in-place.
   * @details If `true`, then `x` will be modified in place if
   *   it is normalized. In case when `x` is dense, it will be both
   *   centered and scaled. If `x` is sparse, it will be only scaled.
   * @param modify_x Whether to modfiy `x` in place or not
   */
  void setModifyX(const bool modify_x);

  /**
   * @brief Sets tolerance in deviance change for early stopping.
   * @param dev_change_tol The tolerance for the change in deviance.
   */
  void setDevChangeTol(const double dev_change_tol);

  /**
   * @brief Sets tolerance in deviance change for early stopping.
   * @param dev_ratio_tol The tolerance for the dev ratio. If the deviance
   * exceeds this value, the path will terminate.
   */
  void setDevRatioTol(const double dev_ratio_tol);

  /**
   * @brief Sets tolerance in deviance change for early stopping.
   * @param max_clusters The maximum number of clusters. SLOPE
   * can (theoretically) select at most select min(n, p) clusters (unique
   * non-zero betas). By default, this is set to -1, which means that the number
   * of clusters will be automatically set to the number of observations + 1.
   */
  void setMaxClusters(const int max_clusters);

  /**
   * @brief Sets the center points for feature normalization.
   * @param type Type of centering, one of: "mean", "none"
   */
  void setCentering(const std::string& type);

  /**
   * @brief Sets the center points for feature normalization.
   * @param x_centers Vector containing center values for each feature
   * Used in feature normalization: x_normalized = (x - center) / scale
   */
  void setCentering(const Eigen::VectorXd& x_centers);

  /**
   * @brief Sets the scaling type
   * @param type Type of scaling, one of: "sd", "l1", "none"
   */
  void setScaling(const std::string& type);

  /**
   * @brief Toggles collection of diagnostics.
   * @param collect_diagnostics Whether to collect diagnostics, i.e.
   * dual gap, objective value, etc. These will be stored in the SlopeFit and
   * SlopePath objects.
   */
  void setDiagnostics(const bool collect_diagnostics);

  /**
   * @brief Sets the scaling factors for feature normalization.
   * @param x_scales Vector containing scale values for each feature
   * Used in feature normalization: x_normalized = (x - center) / scale
   */
  void setScaling(const Eigen::VectorXd& x_scales);

  /**
   * @brief Get currently defined loss type
   * @return The loss type
   */
  const std::string& getLossType();

  /**
   * @brief Computes SLOPE regression solution path for multiple alpha and
   * lambda values
   *
   * @tparam T Matrix type for feature input (supports dense or sparse matrices)
   * @param x Feature matrix of size n x p
   * @param y_in Response matrix of size n x m
   * @param alpha Sequence of mixing parameters for elastic net regularization
   * @param lambda Sequence of regularization parameters (if empty, computed
   * automatically)
   * @return SlopePath object containing full solution path and optimization
   * metrics
   *
   * Fits SLOPE models for each combination of alpha and lambda values, storing
   * all solutions and optimization metrics in a SlopePath object.
   */
  template<typename T>
  SlopePath path(T& x,
                 const Eigen::MatrixXd& y_in,
                 Eigen::ArrayXd alpha = Eigen::ArrayXd::Zero(0),
                 Eigen::ArrayXd lambda = Eigen::ArrayXd::Zero(0));

  /**
   * @brief Fits a single SLOPE regression model for given alpha and lambda
   * values
   *
   * @tparam T Matrix type for feature input (supports dense or sparse matrices)
   * @param x Feature matrix of size n x p
   * @param y_in Response matrix of size n x m
   * @param alpha Mixing parameter for elastic net regularization
   * @param lambda Vector of regularization parameters (if empty, computed
   * automatically)
   * @return SlopeFit Object containing fitted model and optimization metrics
   *
   * Fits a single SLOPE model with specified regularization parameters,
   * returning coefficients and optimization details in a SlopeFit object.
   */
  template<typename T>
  SlopeFit fit(T& x,
               const Eigen::MatrixXd& y_in,
               const double alpha = 1.0,
               Eigen::ArrayXd lambda = Eigen::ArrayXd::Zero(0));

private:
  // Parameters
  bool intercept;
  bool modify_x;
  bool update_clusters;
  bool collect_diagnostics;
  bool return_clusters;
  double alpha_min_ratio;
  double dev_change_tol;
  double dev_ratio_tol;
  double learning_rate_decr;
  double q;
  double theta1;
  double theta2;
  double tol;
  int max_it;
  int path_length;
  int cd_iterations;
  std::optional<int> max_clusters;
  std::string lambda_type;
  std::string centering_type;
  std::string scaling_type;
  std::string loss_type;
  std::string screening_type;
  std::string solver_type;

  // Data
  Eigen::VectorXd x_centers;
  Eigen::VectorXd x_scales;
};

} // namespace slope
