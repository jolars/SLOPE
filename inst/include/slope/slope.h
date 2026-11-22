/**
 * @file
 * @brief SLOPE (Sorted L-One Penalized Estimation) optimization
 */

#pragma once

#include "clusters.h"
#include "constants.h"
#include "diagnostics.h"
#include "estimate_alpha.h"
#include "logger.h"
#include "losses/loss.h"
#include "losses/setup_loss.h"
#include "math.h"
#include "normalize.h"
#include "observation_weights.h"
#include "regularization_sequence.h"
#include "screening.h"
#include "slope_fit.h"
#include "slope_path.h"
#include "solvers/hybrid_cd.h"
#include "solvers/setup_solver.h"
#include "sorted_l1_norm.h"
#include "timer.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>

/** @namespace slope
 *  @brief Namespace containing SLOPE regression implementation
 */
namespace slope {

namespace detail {
struct SlopeAccess;
}

/**
 * @brief Default no-op interrupt checker
 * @return Always returns false
 */
inline bool
defaultInterruptChecker()
{
  return false;
}

/**
 * @brief The SLOPE model.
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
  Slope() = default;

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
   * Disabling cluster updates is intended for benchmarking and experimental
   * use. A coordinate update can then invalidate the ordering and uniqueness
   * assumptions required by later cluster-coordinate updates.
   *
   * @param update_clusters Whether coordinate descent maintains cluster
   * ordering and membership after each update.
   */
  void setUpdateClusters(bool update_clusters);

  /**
   * @brief Sets the return clusters flag.
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
   * @brief Sets the alpha type
   * @details If `"path"`, values of `"alpha"` are automatically selected to
   * span a range from completele penalization to the almost-unpenalized
   * case: the SLOPE path. If `"estimate"`, which is only available for
   * the quadratic loss, the alpha value is estimated by using the
   * scaled L2 norm of the residuals vector from an OLS fit. In case
   * the number of features exceeds the number of observations plus 30,
   * then this is done iteratively starting at the null (or intercept-only)
   * models, and proceeding by refitting slope with these estimations
   * until the set of selected features no longer changes.
   *
   * @param alpha_type The type of alpha, one of `"path"` and `"estimate"`.
   */
  void setAlphaType(const std::string& alpha_type);

  /**
   * @brief Sets the scale of user-supplied alpha values.
   *
   * @param alpha_scale One of "absolute" (the default) or "relative".
   * Relative values passed to fit() and path() are multiplied by alpha_max,
   * the smallest alpha giving zero coefficients for the current data,
   * normalization, observation weights, and lambda sequence. Values must be
   * finite and nonnegative; values above one are allowed. Automatic paths and
   * noise-based alpha estimates are unaffected. Returned fits always store
   * absolute alpha values.
   */
  void setAlphaScale(const std::string& alpha_scale);

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
   * @brief Sets the tolerance value for the relaxed SLOPE solver.
   *
   * @param tol The value to set for the tolerance value. Must be positive.
   */
  void setRelaxTol(double tol);

  /**
   * @brief Sets the maximum number of outer (IRLS) iterations for the relaxed
   * solver.
   *
   * @param max_it The value to set for the maximum number of iterations. Must
   * be positive.
   */
  void setRelaxMaxOuterIterations(int max_it);

  /**
   * @brief Sets the maximum number of inner iterations for the relaxed solver.
   *
   * @param max_it The value to set for the maximum number of iterations. Must
   * be positive.
   */
  void setRelaxMaxInnerIterations(int max_it);

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
   * @brief Sets the frequence of proximal gradient descent steps.
   *
   * @param cd_type The type of coordinate descent to use in the hybrid solver.
   */
  void setHybridCdType(const std::string& cd_type);

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
   *   - "strong": Strong screening rule
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
   * @brief Sets the maximum number of clusters.
   * @param max_clusters The maximum number of clusters. SLOPE
   * can (theoretically) select at most min(n, p) clusters (unique
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
   * @brief Sets the maximum number of iterations for the
   * alpha estimation procedure
   * @details The current stopping criterion for this solver is that the
   * sets of selected features between two consecutive iterations is
   * identical, but users have experienced instances where these sets
   * will cycle back on forth between two sets, which would make the
   * algorithm fail to converge.
   * @param alpha_est_maxit The maximum number of allowed iterations
   */
  void setAlphaEstimationMaxIterations(const int alpha_est_maxit);

  /**
   * @brief Sets the random seed
   *
   * @param seed The value to set for the path length.
   */
  void setRandomSeed(const int seed);

  /**
   * @brief Sets the random seed
   *
   * @param seed The value to set for the path length.
   */
  void setRandomSeed(std::optional<int> seed);

  /**
   * @brief Checks if a random seed is set
   */
  bool hasRandomSeed() const;

  /**
   * @brief Gets the random seed
   */
  int getRandomSeed() const;

  /**
   * @brief Gets the maximum number of iterations allowed for the
   * alpha estimation procedure
   */
  int getAlphaEstimationMaxIterations() const;

  /**
   * @brief Returns the intercept flag
   */
  bool getFitIntercept() const;

  /**
   * @brief Get currently defined loss type
   * @return The loss type as a string
   */
  const std::string& getLossType();

  /**
   * @brief Computes SLOPE regression solution path for multiple alpha and
   * lambda values
   *
   * @tparam T Matrix type for feature input (supports dense or sparse matrices)
   * @param x Feature matrix of size n x p
   * @param y_in Response matrix of size n x m
   * @param alpha Nonincreasing sequence of penalty scales. Empty input
   * generates a path automatically. Values use the scale set by
   * setAlphaScale().
   * @param lambda Sequence of regularization parameters (if empty, computed
   * automatically)
   * @param check_interrupt Optional lambda to check for user interrupt. It runs
   *   periodically during the path fitting.
   * @param observation_weights Optional finite, nonnegative row weights with
   * positive total weight. Empty weights give each row equal weight. Common
   * rescaling leaves the fitted model unchanged.
   * @return SlopePath object containing full solution path and optimization
   * metrics
   *
   * Fits SLOPE models for each combination of alpha and lambda values, storing
   * all solutions and optimization metrics in a SlopePath object.
   */
  template<typename T>
  SlopePath path(
    Eigen::EigenBase<T>& x,
    const Eigen::MatrixXd& y_in,
    Eigen::ArrayXd alpha = Eigen::ArrayXd::Zero(0),
    Eigen::ArrayXd lambda = Eigen::ArrayXd::Zero(0),
    std::function<bool()> check_interrupt = defaultInterruptChecker,
    const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
  {
    return pathImpl(x,
                    y_in,
                    std::move(alpha),
                    std::move(lambda),
                    check_interrupt,
                    observation_weights,
                    false);
  }

private:
  friend struct detail::SlopeAccess;

  template<typename T>
  SlopePath pathImpl(Eigen::EigenBase<T>& x,
                     const Eigen::MatrixXd& y_in,
                     Eigen::ArrayXd alpha,
                     Eigen::ArrayXd lambda,
                     std::function<bool()> check_interrupt,
                     const Eigen::VectorXd& observation_weights,
                     bool response_prepared)
  {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int n = x.rows();
    const int p = x.cols();
    const VectorXd weights =
      detail::normalizeObservationWeights(observation_weights, n);
    const int n_positive = detail::positiveWeightCount(weights);
    if (!detail::uniformObservationWeights(weights) &&
        (this->alpha_type == "estimate" ||
         (lambda.size() == 0 && this->lambda_type == "gaussian"))) {
      throw std::invalid_argument(
        "Nonuniform observation weights are not supported for noise-based "
        "alpha estimation or automatic Gaussian lambda correction");
    }

    const int INTERRUPT_FREQ = 100;
    bool interrupt = false;

    if (n != y_in.rows()) {
      throw std::invalid_argument(
        "x and y_in must have the same number of rows");
    }

    if (!isFinite(x.derived())) {
      throw std::invalid_argument("x must not contain NA, NaN, or Inf values");
    }

    if (!y_in.array().isFinite().all()) {
      throw std::invalid_argument("y must not contain NA, NaN, or Inf values");
    }

    auto jit_normalization = normalize(x.derived(),
                                       this->x_centers,
                                       this->x_scales,
                                       this->centering_type,
                                       this->scaling_type,
                                       this->modify_x,
                                       weights);

    std::unique_ptr<Loss> loss = setupLoss(this->loss_type);

    MatrixXd y = response_prepared ? y_in : loss->preprocessResponse(y_in);

    const int m = y.cols();

    VectorXd beta0 = VectorXd::Zero(m);
    VectorXd beta = VectorXd::Zero(p * m);

    MatrixXd eta = MatrixXd::Zero(n, m); // linear predictor

    if (this->intercept) {
      beta0 = loss->link(detail::weightedMean(y, weights)).transpose();
      eta.rowwise() = beta0.transpose();
    }

    MatrixXd residual = loss->residual(eta, y);
    detail::maskZeroWeightRows(residual, weights);
    VectorXd gradient(beta.size());

    // Path data
    bool user_alpha = alpha.size() > 0;
    bool user_lambda = lambda.size() > 0;

    if (!user_lambda) {
      lambda = lambdaSequence(p * m,
                              this->q,
                              this->lambda_type,
                              n_positive,
                              this->theta1,
                              this->theta2);
    } else {
      if (lambda.size() != beta.size()) {
        throw std::invalid_argument(
          "lambda must be the same length as the number of coefficients");
      }
      if (lambda.minCoeff() < 0) {
        throw std::invalid_argument("lambda must be non-negative");
      }
      if (!lambda.isFinite().all()) {
        throw std::invalid_argument("lambda must be finite");
      }
      // Check that lambda is in decreasing order
      for (int i = 1; i < lambda.size(); ++i) {
        if (lambda(i) > lambda(i - 1)) {
          throw std::invalid_argument("lambda must be in decreasing order");
        }
      }
    }

    // Setup the regularization sequence and path
    SortedL1Norm sl1_norm;

    // TODO: Make this part of the slope class
    auto solver = setupSolver(this->solver_type,
                              this->loss_type,
                              jit_normalization,
                              this->intercept,
                              this->update_clusters,
                              this->cd_iterations,
                              this->cd_type,
                              this->random_seed);

    updateGradient(gradient,
                   x.derived(),
                   residual,
                   this->x_centers,
                   this->x_scales,
                   weights,
                   jit_normalization);

    int alpha_max_ind = whichMax(gradient.cwiseAbs());
    double alpha_max =
      lambda.maxCoeff() == 0.0 ? 0.0 : sl1_norm.dualNorm(gradient, lambda);
    const double numerical_stationarity_tol =
      std::sqrt(std::numeric_limits<double>::epsilon()) *
      std::max(1.0, gradient.cwiseAbs().maxCoeff());

    if (alpha_type == "path" ||
        (alpha_type == "estimate" && alpha_estimate != 1)) {
      const double min_ratio = alpha_min_ratio < 0
                                 ? (n_positive > gradient.size() ? 1e-4 : 1e-2)
                                 : alpha_min_ratio;
      alpha = regularizationPath(alpha, path_length, min_ratio, alpha_max);
      if (user_alpha && alpha_scale == "relative") {
        alpha *= alpha_max;
        if (!alpha.isFinite().all()) {
          throw std::invalid_argument(
            "Relative alpha values must produce finite absolute penalties");
        }
      }
      path_length = alpha.size();
    } else if (alpha_type == "estimate" && alpha_estimate == -1) {
      if (loss_type != "quadratic") {
        throw std::invalid_argument("Automatic alpha estimation is only "
                                    "available for the quadratic loss");
      }
    }

    // Screening setup
    std::unique_ptr<ScreeningRule> screening_rule =
      createScreeningRule(this->screening_type);
    std::vector<int> working_set =
      screening_rule->initialize(static_cast<int>(beta.size()), alpha_max_ind);

    // Path variables
    double null_deviance = loss->deviance(eta, y, weights);
    double dev_prev = null_deviance;

    Timer timer;

    double alpha_prev = std::max(alpha_max, alpha(0));

    std::vector<SlopeFit> fits;
    std::optional<MatrixXd> quadratic_ols_dual_point;

    // Regularization path loop
    for (int path_step = 0; path_step < this->path_length; ++path_step) {
      // Check for interrupt at the start of each path step
      bool local_interrupt = false;
#ifdef _OPENMP
#pragma omp critical(check_interrupt)
#endif
      {
        local_interrupt = check_interrupt();
      }
      if (local_interrupt) {
        interrupt = true;
        break;
      }

      double alpha_curr = alpha(path_step);

      assert(alpha_curr <= alpha_prev && "Alpha must be decreasing");

      Eigen::ArrayXd lambda_curr = alpha_curr * lambda;
      Eigen::ArrayXd lambda_prev = alpha_prev * lambda;

      const bool near_unregularized =
        lambda_curr.maxCoeff() == 0.0 ||
        (alpha_max > 0.0 &&
         alpha_curr <=
           std::sqrt(std::numeric_limits<double>::epsilon()) * alpha_max);
      // The QR residual avoids amplifying roundoff in the loss gradient when
      // the penalty is below numerical resolution. It remains only a candidate
      // until computeDualFromPoint checks and, if necessary, scales it below.
      if (loss_type == "quadratic" && near_unregularized &&
          !quadratic_ols_dual_point.has_value()) {
        quadratic_ols_dual_point =
          detail::quadraticOlsDualPoint(x.derived(),
                                        y,
                                        this->x_centers,
                                        this->x_scales,
                                        jit_normalization,
                                        this->intercept,
                                        weights);
      }

      std::optional<double> quadratic_ols_dual;
      if (quadratic_ols_dual_point.has_value()) {
        quadratic_ols_dual = lambda_curr.maxCoeff() == 0.0
                               ? 0.0
                               : computeDualFromPoint(beta,
                                                      *quadratic_ols_dual_point,
                                                      loss,
                                                      sl1_norm,
                                                      lambda_curr,
                                                      x.derived(),
                                                      y,
                                                      this->x_centers,
                                                      this->x_scales,
                                                      jit_normalization,
                                                      weights);
      }

      std::vector<double> duals, primals, time;
      timer.start();

      // Update gradient for the full set
      // TODO: Only update for non-working set since gradient is updated before
      // the convergence check in the inner loop for the working set
      updateGradient(gradient,
                     x.derived(),
                     residual,
                     x_centers,
                     x_scales,
                     weights,
                     jit_normalization);

      screening_rule->screen(
        working_set, gradient, lambda_curr, lambda_prev, beta);

      int working_set_it = 0;
      int total_it = 0;
      for (; total_it < this->max_it; ++working_set_it, ++total_it) {
        // Compute primal, dual, and gap
        residual = loss->residual(eta, y);
        detail::maskZeroWeightRows(residual, weights);
        updateGradient(gradient,
                       x.derived(),
                       residual,
                       working_set,
                       this->x_centers,
                       this->x_scales,
                       weights,
                       jit_normalization);

        double primal = loss->loss(eta, y, weights) +
                        sl1_norm.eval(beta(working_set),
                                      lambda_curr.head(working_set.size()));

        // The OLS bound can be loose on correlated designs, so it may certify
        // convergence only after the loss itself is numerically stationary.
        const bool numerically_stationary =
          gradient(working_set).cwiseAbs().maxCoeff() <=
            numerical_stationarity_tol &&
          (!this->intercept ||
           detail::weightedMean(residual, weights).cwiseAbs().maxCoeff() <=
             numerical_stationarity_tol);

        MatrixXd dual_point = loss->dualPoint(eta, y, this->intercept, weights);
        MatrixXd theta = dual_point;
        VectorXd dual_gradient = VectorXd::Zero(beta.size());
        updateGradient(dual_gradient,
                       x.derived(),
                       theta,
                       working_set,
                       this->x_centers,
                       this->x_scales,
                       weights,
                       jit_normalization);
        double dual_norm =
          sl1_norm.dualNorm(dual_gradient(working_set),
                            lambda_curr.head(working_set.size()),
                            constants::MAX_DIV);
        theta.array() /= std::max(1.0, dual_norm);
        double dual = loss->dual(theta, y, weights);
        if (quadratic_ols_dual.has_value() && numerically_stationary) {
          dual = std::max(dual, *quadratic_ols_dual);
        }

        double full_dual = std::numeric_limits<double>::quiet_NaN();

        if (collect_diagnostics) {
          timer.pause();
          full_dual = computeDualFromPoint(beta,
                                           dual_point,
                                           loss,
                                           sl1_norm,
                                           lambda_curr,
                                           x.derived(),
                                           y,
                                           this->x_centers,
                                           this->x_scales,
                                           jit_normalization,
                                           weights);
          if (quadratic_ols_dual.has_value() && numerically_stationary) {
            full_dual = std::max(full_dual, *quadratic_ols_dual);
          }
          timer.resume();

          time.emplace_back(timer.elapsed());
          primals.emplace_back(primal);
          duals.emplace_back(full_dual);
        }

        double dual_gap = primal - dual;

        assert(dual_gap > -1e-6 && "Dual gap should be positive");

        double tol_scaled = (std::abs(primal) + constants::EPSILON) * this->tol;
        const bool unregularized_stationary = loss_type == "quadratic" &&
                                              lambda_curr.maxCoeff() == 0.0 &&
                                              numerically_stationary;

        if (dual_gap <= tol_scaled || unregularized_stationary ||
            working_set_it + 1 == this->max_it) {
          MatrixXd screening_residual = residual;
          for (int i = 0; i < n; ++i) {
            if (weights(i) == 0.0) {
              screening_residual.row(i).setZero();
            } else {
              screening_residual.row(i) *= weights(i);
            }
          }
          bool no_violations =
            screening_rule->checkKktViolations(gradient,
                                               beta,
                                               lambda_curr,
                                               working_set,
                                               x.derived(),
                                               screening_residual,
                                               this->x_centers,
                                               this->x_scales,
                                               jit_normalization);
          if (no_violations) {
            if (!std::isfinite(full_dual)) {
              full_dual = computeDualFromPoint(beta,
                                               dual_point,
                                               loss,
                                               sl1_norm,
                                               lambda_curr,
                                               x.derived(),
                                               y,
                                               this->x_centers,
                                               this->x_scales,
                                               jit_normalization,
                                               weights);
              if (quadratic_ols_dual.has_value() && numerically_stationary) {
                full_dual = std::max(full_dual, *quadratic_ols_dual);
              }
            }
            if (primal - full_dual <= tol_scaled || unregularized_stationary) {
              break;
            }
          } else {
            working_set_it = 0; // Restart if there are KKT violations
          }
        }

        if (working_set_it % INTERRUPT_FREQ == 0) {
          bool local_interrupt = false;
#ifdef _OPENMP
#pragma omp critical(check_interrupt)
#endif
          {
            local_interrupt = check_interrupt();
          }
          if (local_interrupt) {
            interrupt = true;
            break;
          }
        }

        solver->run(beta0,
                    beta,
                    eta,
                    lambda_curr,
                    loss,
                    sl1_norm,
                    gradient,
                    working_set,
                    x.derived(),
                    this->x_centers,
                    this->x_scales,
                    y,
                    weights);
      }

      if (total_it == this->max_it) {
        WarningLogger::addWarning(
          WarningCode::MAXIT_REACHED,
          "Maximum number of iterations reached at step = " +
            std::to_string(path_step) + ".");
      }

      alpha_prev = alpha_curr;

      // Compute early stopping criteria
      double dev = loss->deviance(eta, y, weights);
      double dev_ratio = 1 - dev / null_deviance;
      double dev_change = path_step == 0 ? 1.0 : 1 - dev / dev_prev;
      dev_prev = dev;

      Clusters clusters;

      if (return_clusters) {
        clusters.update(beta);
      }

      SlopeFit fit{ beta0,
                    beta.reshaped(p, m).sparseView(),
                    clusters,
                    alpha_curr,
                    lambda,
                    dev,
                    null_deviance,
                    primals,
                    duals,
                    time,
                    total_it,
                    this->centering_type,
                    this->scaling_type,
                    this->intercept,
                    this->x_centers,
                    this->x_scales };

      fits.emplace_back(std::move(fit));

      if (interrupt) {
        break;
      }

      if (!user_alpha) {
        int n_unique = unique(beta.cwiseAbs()).size();
        if (dev_ratio > dev_ratio_tol || dev_change < dev_change_tol ||
            n_unique >= this->max_clusters.value_or(n_positive + 1)) {
          break;
        }
      }
    }

    return fits;
  }

public:
  /**
   * @brief Fits a single SLOPE regression model for given alpha and lambda
   * values
   *
   * @tparam T Matrix type for feature input (supports dense or sparse matrices)
   * @param x Feature matrix of size n x p
   * @param y_in Response matrix of size n x m
   * @param alpha Nonnegative penalty scale, interpreted according to
   * setAlphaScale(). Defaults to 1.0, giving zero coefficients when the scale
   * is relative and the lambda sequence has a positive entry.
   * @param lambda Vector of regularization parameters (if empty, computed
   * automatically)
   * @param check_interrupt Optional lambda to check for user interrupt. It runs
   *  periodically during the fitting.
   * @return SlopeFit Object containing fitted model and optimization metrics
   *
   * Fits a single SLOPE model with specified regularization parameters,
   * returning coefficients and optimization details in a SlopeFit object.
   */
  template<typename T>
  SlopeFit fit(Eigen::EigenBase<T>& x,
               const Eigen::MatrixXd& y_in,
               const double alpha = 1.0,
               Eigen::ArrayXd lambda = Eigen::ArrayXd::Zero(0),
               std::function<bool()> check_interrupt = defaultInterruptChecker,
               const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
  {
    Eigen::ArrayXd alpha_arr(1);
    alpha_arr(0) = alpha;
    SlopePath res =
      path(x, y_in, alpha_arr, lambda, check_interrupt, observation_weights);

    return { res(0) };
  };

  /**
   * @brief Estimates the regularization parameter alpha for SLOPE regression
   *
   * This function implements an algorithm to estimate an appropriate
   * regularization parameter (alpha) for SLOPE, which is a generalization of
   * the lasso. When n >= p + 30, it directly estimates alpha from OLS
   * residuals. Otherwise, it uses an iterative procedure that alternates
   * between estimating alpha and fitting the SLOPE model.
   *
   * The iterative procedure works by:
   * 1. Starting with an empty set of selected variables
   * 2. Estimating alpha based on the selected variables
   * 3. Fitting a SLOPE model with that alpha
   * 4. Updating the selected variables based on non-zero coefficients
   * 5. Repeating until convergence or maximum iterations reached
   *
   * @tparam MatrixType The type of matrix used to store the design matrix
   * @param x Design matrix with n observations and p predictors
   * @param y Response matrix
   * @param check_interrupt Optional lambda to check for user interrupt
   * @param observation_weights Optional row weights. Only constant positive
   * weights are supported by noise-based alpha estimation.
   * @return A SlopePath object containing the fitted model with estimated alpha
   * @throws std::runtime_error If maximum iterations reached or if too many
   * variables selected
   */
  template<typename T>
  SlopePath estimateAlpha(
    Eigen::EigenBase<T>& x,
    Eigen::MatrixXd& y,
    std::function<bool()> check_interrupt = defaultInterruptChecker,
    const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
  {
    const auto weights =
      detail::normalizeObservationWeights(observation_weights, x.rows());
    if (!detail::uniformObservationWeights(weights)) {
      throw std::invalid_argument("Nonuniform observation weights are not "
                                  "supported for noise-based alpha estimation");
    }
    int n = x.rows();
    int p = x.cols();

    // Create a copy with alpha type set to path to avoid recursion
    Slope model_copy = *this;
    model_copy.setAlphaType("path");
    model_copy.setAlphaScale("absolute");

    std::vector<int> selected;
    Eigen::ArrayXd alpha(1);
    SlopePath result;

    // Estimate the noise level, if possible
    if (n >= p + 30) {
      alpha(0) = estimateNoise(x, y, this->intercept) / n;
      this->alpha_estimate = alpha(0);
      result = model_copy.path(
        x, y, alpha, Eigen::ArrayXd::Zero(0), check_interrupt, weights);
    } else {
      for (int it = 0; it < this->alpha_est_maxit; ++it) {
        T x_selected = subsetCols(x.derived(), selected);

        std::vector<int> selected_prev = selected;
        selected.clear();

        alpha(0) = estimateNoise(x_selected, y, this->intercept) / n;
        this->alpha_estimate = alpha(0);

        result = model_copy.path(
          x, y, alpha, Eigen::ArrayXd::Zero(0), check_interrupt, weights);
        auto coefs = result.getCoefs().back();

        for (typename Eigen::SparseMatrix<double>::InnerIterator it(coefs, 0);
             it;
             ++it) {
          selected.emplace_back(it.row());
        }

        if (selected == selected_prev) {
          return result;
        }

        if (static_cast<int>(selected.size()) >= n + this->intercept) {
          throw std::runtime_error(
            "selected >= n - 1 variables, cannot estimate variance");
        }
      }

      slope::WarningLogger::addWarning(
        slope::WarningCode::MAXIT_REACHED,
        "Maximum iterations reached in alpha estimation");
    }

    return result;
  }

  /**
   * @brief Relaxes a fitted SLOPE model
   *
   * @tparam T Matrix type for feature input (supports dense or sparse matrices)
   * @param fit Previously fitted SLOPE model containing coefficient estimates
   * @param x Original feature matrix of size n x p, before normalization.
   * Relaxation uses just-in-time normalization and does not modify x.
   * @param y_in Response vector of size n
   * @param gamma Relaxation parameter, proportion of SLOPE-penalized fit. Must
   * be between 0 and 1. Default is 0.0 which means fully relaxed.
   * @param beta0 Warm start intercept values (optional)
   * @param beta Warm start coefficient values (optional)
   * @param observation_weights Optional row weights. Supply the fitting weights
   * again when relaxing a weighted fit.
   * @return SlopeFit Object containing the relaxed model with unpenalized
   * coefficients
   */
  template<typename T>
  SlopeFit relax(const SlopeFit& fit,
                 T& x,
                 const Eigen::VectorXd& y_in,
                 const double gamma = 0.0,
                 Eigen::VectorXd beta0 = Eigen::VectorXd(0),
                 Eigen::VectorXd beta = Eigen::VectorXd(0),
                 const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
  {
    return relaxImpl(fit,
                     x,
                     y_in,
                     gamma,
                     std::move(beta0),
                     std::move(beta),
                     observation_weights,
                     false);
  }

private:
  struct RelaxationState
  {
    Eigen::SparseMatrix<int> pattern;
    Eigen::VectorXd beta0;
    Eigen::VectorXd beta;
    bool initialized = false;
  };

  template<typename T>
  SlopeFit relaxImpl(const SlopeFit& fit,
                     T& x,
                     const Eigen::MatrixXd& y_in,
                     double gamma,
                     Eigen::VectorXd beta0,
                     Eigen::VectorXd beta,
                     const Eigen::VectorXd& observation_weights,
                     bool response_prepared,
                     RelaxationState* state = nullptr)
  {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    int n = x.rows();
    int p = x.cols();
    const VectorXd weights =
      detail::normalizeObservationWeights(observation_weights, n);
    if (y_in.rows() != n || !isFinite(x) || !y_in.allFinite()) {
      throw std::invalid_argument(
        "Relaxation requires finite x and y with matching rows");
    }
    if (!std::isfinite(gamma) || gamma < 0.0 || gamma > 1.0) {
      throw std::invalid_argument("gamma must be between zero and one");
    }
    if (gamma == 1.0) {
      return fit;
    }

    if (beta0.size() == 0) {
      beta0 = fit.getIntercepts(false);
    }

    if (beta.size() == 0) {
      beta = MatrixXd(fit.getCoefs(false)).reshaped();
    }

    auto working_set = activeSet(beta);
    bool reuse_solution = false;
    if (state != nullptr) {
      auto pattern = patternMatrix(beta);
      reuse_solution =
        state->initialized && pattern.rows() == state->pattern.rows() &&
        pattern.cols() == state->pattern.cols() &&
        (pattern.cols() == 0 || pattern.isApprox(state->pattern));
      if (reuse_solution) {
        beta0 = state->beta0;
        beta = state->beta;
      } else {
        state->pattern = std::move(pattern);
      }
    }

    double alpha = fit.getAlpha();

    Timer timer;

    std::vector<double> primals, duals, time;
    timer.start();

    auto jit_normalization = normalize(
      x, x_centers, x_scales, centering_type, scaling_type, false, weights);

    bool update_clusters = false;

    std::unique_ptr<Loss> loss = setupLoss(this->loss_type);

    MatrixXd y = response_prepared ? y_in : loss->preprocessResponse(y_in);

    int m = y.cols();

    MatrixXd eta;
    int passes = 0;

    if (!reuse_solution) {
      Eigen::ArrayXd lambda_cumsum_relax = Eigen::ArrayXd::Zero(p * m + 1);

      eta = linearPredictor(x,
                            working_set,
                            beta0,
                            beta,
                            x_centers,
                            x_scales,
                            jit_normalization,
                            intercept);
      VectorXd gradient = VectorXd::Zero(p * m);
      MatrixXd residual(n, m);
      MatrixXd working_residual(n, m);

      MatrixXd w = MatrixXd::Ones(n, m);
      MatrixXd w_ones = weights.replicate(1, m);
      MatrixXd z = y;

      slope::Clusters clusters(beta);

      std::mt19937 rng;

      if (random_seed.has_value()) {
        rng.seed(*random_seed);
      } else {
        rng.seed(std::random_device{}());
      }

      for (int irls_it = 0; irls_it < max_it_outer_relax; irls_it++) {
        residual = loss->residual(eta, y);
        detail::maskZeroWeightRows(residual, weights);

        if (collect_diagnostics) {
          primals.push_back(loss->loss(eta, y, weights));
          duals.push_back(0.0);
          time.push_back(timer.elapsed());
        }

        Eigen::VectorXd cluster_gradient = clusterGradient(beta,
                                                           residual,
                                                           clusters,
                                                           x,
                                                           w_ones,
                                                           x_centers,
                                                           x_scales,
                                                           jit_normalization);

        double norm_grad = cluster_gradient.lpNorm<Eigen::Infinity>();
        if (intercept) {
          norm_grad = std::max(
            norm_grad,
            detail::weightedMean(residual, weights).cwiseAbs().maxCoeff());
        }

        if (norm_grad < tol_relax) {
          break;
        }

        loss->updateWeightsAndWorkingResponse(w, z, eta, y, weights);
        working_residual = eta - z;
        detail::maskZeroWeightRows(working_residual, weights);
        const VectorXd weight_sums = w.colwise().sum().transpose();

        for (int inner_it = 0; inner_it < max_it_inner_relax; ++inner_it) {
          passes++;

          double max_abs_gradient = coordinateDescent(beta0,
                                                      beta,
                                                      working_residual,
                                                      clusters,
                                                      lambda_cumsum_relax,
                                                      x,
                                                      w,
                                                      weight_sums,
                                                      x_centers,
                                                      x_scales,
                                                      intercept,
                                                      jit_normalization,
                                                      update_clusters,
                                                      rng,
                                                      cd_type);

          if (max_abs_gradient < tol_relax) {
            break;
          }
        }

        eta = working_residual + z;

        if (irls_it == max_it_outer_relax) {
          WarningLogger::addWarning(
            WarningCode::MAXIT_REACHED,
            "Maximum number of IRLS iterations reached.");
        }
      }
    }

    if (state != nullptr && !reuse_solution) {
      // Cache before blending so later steps use their own penalized fit.
      state->beta0 = beta0;
      state->beta = beta;
      state->initialized = true;
    }

    if (gamma > 0) {
      Eigen::VectorXd old_coefs = MatrixXd(fit.getCoefs(false)).reshaped();
      Eigen::VectorXd old_intercept = fit.getIntercepts(false);
      beta = (1 - gamma) * beta + gamma * old_coefs;
      beta0 = (1 - gamma) * beta0 + gamma * old_intercept;
    }

    eta = linearPredictor(x,
                          working_set,
                          beta0,
                          beta,
                          x_centers,
                          x_scales,
                          jit_normalization,
                          intercept);
    const double dev = loss->deviance(eta, y, weights);
    MatrixXd null_eta = MatrixXd::Zero(n, m);
    if (intercept) {
      null_eta.rowwise() = loss->link(detail::weightedMean(y, weights)).row(0);
    }
    const double null_dev = loss->deviance(null_eta, y, weights);
    Clusters clusters(beta);

    SlopeFit fit_out{ beta0,           beta.reshaped(p, m).sparseView(),
                      clusters,        alpha,
                      fit.getLambda(), dev,
                      null_dev,        primals,
                      duals,           time,
                      passes,          centering_type,
                      scaling_type,    intercept,
                      x_centers,       x_scales };

    return fit_out;
  }

  template<typename T>
  SlopePath relaxPathImpl(const SlopePath& path,
                          T& x,
                          const Eigen::MatrixXd& y,
                          double gamma,
                          const Eigen::VectorXd& observation_weights,
                          bool response_prepared)
  {
    std::vector<SlopeFit> fits;
    fits.reserve(path.size());
    RelaxationState state;
    for (size_t i = 0; i < path.size(); ++i) {
      fits.emplace_back(relaxImpl(path(i),
                                  x,
                                  y,
                                  gamma,
                                  {},
                                  {},
                                  observation_weights,
                                  response_prepared,
                                  &state));
    }
    return fits;
  }

public:
  /**
   * @brief Relaxes a fitted SLOPE path
   *
   * Reuses the unpenalized fit at consecutive steps with the same signed
   * cluster pattern. Reused steps have zero optimization passes and empty
   * diagnostic histories.
   *
   * @tparam T Matrix type for feature input (supports dense or sparse matrices)
   * @param path Previously fitted SLOPE path
   * @param x Feature matrix of size n x p
   * @param y Response vector of size n
   * @param gamma Relaxation parameter, proportion of SLOPE-penalized fit. Must
   * be between 0 and 1. Default is 0.0 which means fully relaxed.
   * @param observation_weights Optional row weights, as in the fit overload.
   * @return SlopePath Object containing the relaxed model with unpenalized
   * coefficients
   */
  template<typename T>
  SlopePath relax(
    const SlopePath& path,
    T& x,
    const Eigen::VectorXd& y,
    const double gamma = 0.0,
    const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
  {
    return relaxPathImpl(path, x, y, gamma, observation_weights, false);
  }

private:
  // Parameters
  bool collect_diagnostics = false;
  bool intercept = true;
  bool modify_x = false;
  bool return_clusters = true;
  bool update_clusters = true;
  double alpha_min_ratio = -1;
  double dev_change_tol = 1e-5;
  double dev_ratio_tol = 0.999;
  double learning_rate_decr = 0.5;
  double q = 0.1;
  double theta1 = 1.0;
  double theta2 = 0.5;
  double tol = 1e-4;
  double tol_relax = 1e-4;
  double alpha_estimate = -1;
  int alpha_est_maxit = 1000;
  int cd_iterations = 10;
  int max_it = 1e5;
  int max_it_inner_relax = 1e5;
  int max_it_outer_relax = 50;
  int path_length = 100;
  std::optional<int> max_clusters = std::nullopt;
  std::optional<int> random_seed = 0;
  std::string alpha_scale = "absolute";
  std::string alpha_type = "path";
  std::string cd_type = "permuted";
  std::string centering_type = "mean";
  std::string lambda_type = "bh";
  std::string loss_type = "quadratic";
  std::string scaling_type = "sd";
  std::string screening_type = "strong";
  std::string solver_type = "auto";

  // Data
  Eigen::VectorXd x_centers;
  Eigen::VectorXd x_scales;
};

namespace detail {

/** CV prepares responses once to keep class columns stable across folds. */
struct SlopeAccess
{
  template<typename T>
  static SlopePath path(Slope& model,
                        Eigen::EigenBase<T>& x,
                        const Eigen::MatrixXd& y,
                        const Eigen::ArrayXd& alphas,
                        const Eigen::VectorXd& weights)
  {
    return model.pathImpl(
      x, y, alphas, {}, defaultInterruptChecker, weights, true);
  }

  template<typename T>
  static SlopePath relax(Slope& model,
                         const SlopePath& path,
                         T& x,
                         const Eigen::MatrixXd& y,
                         double gamma,
                         const Eigen::VectorXd& weights)
  {
    return model.relaxPathImpl(path, x, y, gamma, weights, true);
  }
};

} // namespace detail

} // namespace slope
