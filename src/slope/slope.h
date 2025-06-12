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
#include <memory>
#include <numeric>
#include <optional>

/** @namespace slope
 *  @brief Namespace containing SLOPE regression implementation
 */
namespace slope {

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
   * @param update_clusters Selects whether the coordinate descent keeps the
   * clusters updated.
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
  SlopePath path(Eigen::EigenBase<T>& x,
                 const Eigen::MatrixXd& y_in,
                 Eigen::ArrayXd alpha = Eigen::ArrayXd::Zero(0),
                 Eigen::ArrayXd lambda = Eigen::ArrayXd::Zero(0))
  {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int n = x.rows();
    const int p = x.cols();

    if (n != y_in.rows()) {
      throw std::invalid_argument(
        "x and y_in must have the same number of rows");
    }

    auto jit_normalization = normalize(x.derived(),
                                       this->x_centers,
                                       this->x_scales,
                                       this->centering_type,
                                       this->scaling_type,
                                       this->modify_x);

    std::unique_ptr<Loss> loss = setupLoss(this->loss_type);

    MatrixXd y = loss->preprocessResponse(y_in);

    const int m = y.cols();

    std::vector<int> full_set(p * m);
    std::iota(full_set.begin(), full_set.end(), 0);

    VectorXd beta0 = VectorXd::Zero(m);
    VectorXd beta = VectorXd::Zero(p * m);

    MatrixXd eta = MatrixXd::Zero(n, m); // linear predictor

    if (this->intercept) {
      beta0 = loss->link(y.colwise().mean()).transpose();
      eta.rowwise() = beta0.transpose();
    }

    MatrixXd residual = loss->residual(eta, y);
    VectorXd gradient(beta.size());

    // Path data
    bool user_alpha = alpha.size() > 0;
    bool user_lambda = lambda.size() > 0;

    if (!user_lambda) {
      lambda = lambdaSequence(
        p * m, this->q, this->lambda_type, n, this->theta1, this->theta2);
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
    }

    // Setup the regularization sequence and path
    SortedL1Norm sl1_norm;

    // TODO: Make this part of the slope class
    auto solver = setupSolver(this->solver_type,
                              this->loss_type,
                              jit_normalization,
                              this->intercept,
                              this->update_clusters,
                              this->cd_iterations);

    updateGradient(gradient,
                   x.derived(),
                   residual,
                   full_set,
                   this->x_centers,
                   this->x_scales,
                   Eigen::VectorXd::Ones(n),
                   jit_normalization);

    int alpha_max_ind = whichMax(gradient.cwiseAbs());
    double alpha_max = sl1_norm.dualNorm(gradient, lambda);

    if (alpha_type == "path" ||
        (alpha_type == "estimate" && alpha_estimate != 1)) {
      if (alpha_min_ratio < 0) {
        alpha_min_ratio = n > gradient.size() ? 1e-4 : 1e-2;
      }

      alpha =
        regularizationPath(alpha, path_length, alpha_min_ratio, alpha_max);
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
      screening_rule->initialize(full_set, alpha_max_ind);

    // Path variables
    double null_deviance = loss->deviance(eta, y);
    double dev_prev = null_deviance;

    Timer timer;

    double alpha_prev = std::max(alpha_max, alpha(0));

    std::vector<SlopeFit> fits;

    // Regularization path loop
    for (int path_step = 0; path_step < this->path_length; ++path_step) {
      double alpha_curr = alpha(path_step);

      assert(alpha_curr <= alpha_prev && "Alpha must be decreasing");

      Eigen::ArrayXd lambda_curr = alpha_curr * lambda;
      Eigen::ArrayXd lambda_prev = alpha_prev * lambda;

      std::vector<double> duals, primals, time;
      timer.start();

      // Update gradient for the full set
      // TODO: Only update for non-working set since gradient is updated before
      // the convergence check in the inner loop for the working set
      updateGradient(gradient,
                     x.derived(),
                     residual,
                     full_set,
                     x_centers,
                     x_scales,
                     Eigen::VectorXd::Ones(x.rows()),
                     jit_normalization);

      working_set = screening_rule->screen(
        gradient, lambda_curr, lambda_prev, beta, full_set);

      int it = 0;
      int total_it = 0;
      for (; it < this->max_it; ++it, ++total_it) {
        // Compute primal, dual, and gap
        residual = loss->residual(eta, y);
        updateGradient(gradient,
                       x.derived(),
                       residual,
                       working_set,
                       this->x_centers,
                       this->x_scales,
                       Eigen::VectorXd::Ones(n),
                       jit_normalization);

        double primal = loss->loss(eta, y) +
                        sl1_norm.eval(beta(working_set),
                                      lambda_curr.head(working_set.size()));

        MatrixXd theta = residual;

        // First compute gradient with potential offset for intercept case
        VectorXd dual_gradient = gradient;

        // TODO: Can we avoid this copy? Maybe revert offset afterwards or,
        // alternatively, solve intercept until convergence and then no longer
        // need the offset at all.
        if (this->intercept) {
          VectorXd theta_mean = theta.colwise().mean();
          theta.rowwise() -= theta_mean.transpose();

          offsetGradient(dual_gradient,
                         x.derived(),
                         theta_mean,
                         working_set,
                         this->x_centers,
                         this->x_scales,
                         jit_normalization);
        }

        // Common scaling operation
        double dual_norm = sl1_norm.dualNorm(
          dual_gradient(working_set), lambda_curr.head(working_set.size()));
        theta.array() /= std::max(1.0, dual_norm);

        double dual = loss->dual(theta, y, Eigen::VectorXd::Ones(n));

        if (collect_diagnostics) {
          timer.pause();
          double true_dual = computeDual(beta,
                                         residual,
                                         loss,
                                         sl1_norm,
                                         lambda_curr,
                                         x.derived(),
                                         y,
                                         this->x_centers,
                                         this->x_scales,
                                         jit_normalization,
                                         this->intercept);
          timer.resume();

          time.emplace_back(timer.elapsed());
          primals.emplace_back(primal);
          duals.emplace_back(true_dual);
        }

        double dual_gap = primal - dual;

        assert(dual_gap > -1e-6 && "Dual gap should be positive");

        double tol_scaled = (std::abs(primal) + constants::EPSILON) * this->tol;

        if (dual_gap <= tol_scaled || it == this->max_it) {
          bool no_violations =
            screening_rule->checkKktViolations(gradient,
                                               beta,
                                               lambda_curr,
                                               working_set,
                                               x.derived(),
                                               residual,
                                               this->x_centers,
                                               this->x_scales,
                                               jit_normalization,
                                               full_set);
          if (no_violations) {
            break;
          } else {
            it = 0; // Restart if there are KKT violations
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
                    y);
      }

      if (it == this->max_it) {
        WarningLogger::addWarning(
          WarningCode::MAXIT_REACHED,
          "Maximum number of iterations reached at step = " +
            std::to_string(path_step) + ".");
      }

      alpha_prev = alpha_curr;

      // Compute early stopping criteria
      double dev = loss->deviance(eta, y);
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

      if (!user_alpha) {
        int n_unique = unique(beta.cwiseAbs()).size();
        if (dev_ratio > dev_ratio_tol || dev_change < dev_change_tol ||
            n_unique >= this->max_clusters.value_or(n + 1)) {
          break;
        }
      }
    }

    return fits;
  }

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
  SlopeFit fit(Eigen::EigenBase<T>& x,
               const Eigen::MatrixXd& y_in,
               const double alpha = 1.0,
               Eigen::ArrayXd lambda = Eigen::ArrayXd::Zero(0))
  {
    Eigen::ArrayXd alpha_arr(1);
    alpha_arr(0) = alpha;
    SlopePath res = path(x, y_in, alpha_arr, lambda);

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
   * @return A SlopePath object containing the fitted model with estimated alpha
   * @throws std::runtime_error If maximum iterations reached or if too many
   * variables selected
   */
  template<typename T>
  SlopePath estimateAlpha(Eigen::EigenBase<T>& x, Eigen::MatrixXd& y)
  {
    int n = x.rows();
    int p = x.cols();

    // Create a copy with alpha type set to path to avoid recursion
    Slope model_copy = *this;
    model_copy.setAlphaType("path");

    std::vector<int> selected;
    Eigen::ArrayXd alpha(1);
    SlopePath result;

    // Estimate the noise level, if possible
    if (n >= p + 30) {
      alpha(0) = estimateNoise(x, y, this->intercept) / n;
      this->alpha_estimate = alpha(0);
      result = model_copy.path(x, y, alpha);
    } else {
      for (int it = 0; it < this->alpha_est_maxit; ++it) {
        T x_selected = subsetCols(x.derived(), selected);

        std::vector<int> selected_prev = selected;
        selected.clear();

        alpha(0) = estimateNoise(x_selected, y, this->intercept) / n;
        this->alpha_estimate = alpha(0);

        result = model_copy.path(x, y, alpha);
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
   * @param x Feature matrix of size n x p
   * @param y_in Response vector of size n
   * @param gamma Relaxation parameter, proportion of SLOPE-penalized fit. Must
   * be between 0 and 1. Default is 0.0 which means fully relaxed.
   * @param beta0 Warm start intercept values (optional)
   * @param beta Warm start coefficient values (optional)
   * @return SlopeFit Object containing the relaxed model with unpenalized
   * coefficients
   */
  template<typename T>
  SlopeFit relax(const SlopeFit& fit,
                 T& x,
                 const Eigen::VectorXd& y_in,
                 const double gamma = 0.0,
                 Eigen::VectorXd beta0 = Eigen::VectorXd(0),
                 Eigen::VectorXd beta = Eigen::VectorXd(0))
  {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    int n = x.rows();
    int p = x.cols();

    if (beta0.size() == 0) {
      beta0 = fit.getIntercepts(false);
    }

    if (beta.size() == 0) {
      beta = fit.getCoefs(false);
    }

    double alpha = fit.getAlpha();

    Timer timer;

    std::vector<double> primals, duals, time;
    timer.start();

    auto jit_normalization =
      normalize(x, x_centers, x_scales, centering_type, scaling_type, modify_x);

    bool update_clusters = false;

    std::unique_ptr<Loss> loss = setupLoss(this->loss_type);

    MatrixXd y = loss->preprocessResponse(y_in);

    int m = y.cols();

    Eigen::ArrayXd lambda_relax = Eigen::ArrayXd::Zero(p * m);

    auto working_set = activeSet(beta);

    Eigen::MatrixXd eta = linearPredictor(x,
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
    MatrixXd w_ones = MatrixXd::Ones(n, m);
    MatrixXd z = y;

    Clusters clusters = fit.getClusters();

    int passes = 0;

    for (int irls_it = 0; irls_it < max_it_outer_relax; irls_it++) {
      residual = loss->residual(eta, y);

      if (collect_diagnostics) {
        primals.push_back(loss->loss(eta, y));
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

      if (norm_grad < tol_relax) {
        break;
      }

      loss->updateWeightsAndWorkingResponse(w, z, eta, y);
      working_residual = eta - z;

      for (int inner_it = 0; inner_it < max_it_inner_relax; ++inner_it) {
        passes++;

        double max_abs_gradient = coordinateDescent(beta0,
                                                    beta,
                                                    working_residual,
                                                    clusters,
                                                    lambda_relax,
                                                    x,
                                                    w,
                                                    x_centers,
                                                    x_scales,
                                                    intercept,
                                                    jit_normalization,
                                                    update_clusters);

        if (max_abs_gradient < tol_relax) {
          break;
        }
      }

      eta = working_residual + z;

      if (irls_it == max_it_outer_relax) {
        WarningLogger::addWarning(WarningCode::MAXIT_REACHED,
                                  "Maximum number of IRLS iterations reached.");
      }
    }

    double dev = loss->deviance(eta, y);

    if (gamma > 0) {
      Eigen::VectorXd old_coefs = fit.getCoefs(false);
      Eigen::VectorXd old_intercept = fit.getIntercepts(false);
      beta = (1 - gamma) * beta + gamma * old_coefs;
    }

    SlopeFit fit_out{ beta0,
                      beta.reshaped(p, m).sparseView(),
                      clusters,
                      alpha,
                      lambda_relax,
                      dev,
                      fit.getNullDeviance(),
                      primals,
                      duals,
                      time,
                      passes,
                      centering_type,
                      scaling_type,
                      intercept,
                      x_centers,
                      x_scales };

    return fit_out;
  }

  /**
   * @brief Relaxes a fitted SLOPE path
   *
   * @tparam T Matrix type for feature input (supports dense or sparse matrices)
   * @param path Previously fitted SLOPE path
   * @param x Feature matrix of size n x p
   * @param y Response vector of size n
   * @param gamma Relaxation parameter, proportion of SLOPE-penalized fit. Must
   * be between 0 and 1. Default is 0.0 which means fully relaxed.
   * @return SlopePath Object containing the relaxed model with unpenalized
   * coefficients
   */
  template<typename T>
  SlopePath relax(const SlopePath& path,
                  T& x,
                  const Eigen::VectorXd& y,
                  const double gamma = 0.0)
  {
    std::vector<SlopeFit> fits;

    Eigen::VectorXd beta0 = path(0).getIntercepts(false);
    Eigen::VectorXd beta = path(0).getCoefs(false);

    for (size_t i = 0; i < path.size(); i++) {
      auto relaxed_fit = relax(path(i), x, y, gamma, beta0, beta);

      fits.emplace_back(relaxed_fit);

      // Update warm starts
      // TODO: Maybe be more clever about whether to use the
      // previous values or the regularized estimates and warm starts.
      // Maybe just pick the solution with larger coefficients?
      beta0 = relaxed_fit.getIntercepts(false);
      beta = relaxed_fit.getCoefs(false);
    }

    return fits;
  }

private:
  // Parameters
  bool collect_diagnostics = false;
  bool intercept = true;
  bool modify_x = false;
  bool return_clusters = true;
  bool update_clusters = false;
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
  std::string alpha_type = "path";
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

} // namespace slope
