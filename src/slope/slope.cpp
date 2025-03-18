#include "slope.h"
#include "clusters.h"
#include "constants.h"
#include "diagnostics.h"
#include "estimate_alpha.h"
#include "kkt_check.h"
#include "logger.h"
#include "losses/loss.h"
#include "losses/setup_loss.h"
#include "math.h"
#include "normalize.h"
#include "regularization_sequence.h"
#include "screening.h"
#include "solvers/setup_solver.h"
#include "sorted_l1_norm.h"
#include "timer.h"
#include "utils.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>
#include <numeric>
#include <set>

namespace slope {

template<typename T>
SlopePath
Slope::path(T& x,
            const Eigen::MatrixXd& y_in,
            Eigen::ArrayXd alpha,
            Eigen::ArrayXd lambda)
{
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  const int n = x.rows();
  const int p = x.cols();

  if (n != y_in.rows()) {
    throw std::invalid_argument("x and y_in must have the same number of rows");
  }

  auto jit_normalization = normalize(x,
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
                 x,
                 residual,
                 full_set,
                 this->x_centers,
                 this->x_scales,
                 Eigen::VectorXd::Ones(n),
                 jit_normalization);

  int alpha_max_ind = whichMax(gradient.cwiseAbs());
  double alpha_max = sl1_norm.dualNorm(gradient, lambda);

  if (alpha_type == "path") {
    if (alpha_min_ratio < 0) {
      alpha_min_ratio = n > gradient.size() ? 1e-4 : 1e-2;
    }

    alpha = regularizationPath(alpha, path_length, alpha_min_ratio, alpha_max);
    path_length = alpha.size();
  } else {
    if (loss_type != "quadratic") {
      throw std::invalid_argument(
        "Automatic alpha estimation is only available for the quadratic loss");
    }
    Slope model = *this;
    return estimateAlpha(x, y, model);
  }

  // Screening stuff
  std::vector<int> strong_set, previous_set, working_set, inactive_set;
  if (this->screening_type == "none") {
    working_set = full_set;
  } else {
    working_set = { alpha_max_ind };
  }

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

    if (screening_type == "strong") {
      // TODO: Only update for inactive set, making sure gradient
      // is updated for the active set
      updateGradient(gradient,
                     x,
                     residual,
                     full_set,
                     this->x_centers,
                     this->x_scales,
                     Eigen::VectorXd::Ones(n),
                     jit_normalization);

      previous_set = activeSet(beta);
      strong_set = strongSet(gradient, lambda_curr, lambda_prev);
      strong_set = setUnion(strong_set, previous_set);
      working_set = setUnion(previous_set, { alpha_max_ind });
    }

    int it = 0;
    for (; it < this->max_it; ++it) {
      // Compute primal, dual, and gap
      residual = loss->residual(eta, y);
      updateGradient(gradient,
                     x,
                     residual,
                     working_set,
                     this->x_centers,
                     this->x_scales,
                     Eigen::VectorXd::Ones(n),
                     jit_normalization);

      double primal =
        loss->loss(eta, y) +
        sl1_norm.eval(beta(working_set), lambda_curr.head(working_set.size()));

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
                       x,
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
                                       x,
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

      if (dual_gap <= tol_scaled) {
        if (screening_type == "strong") {
          updateGradient(gradient,
                         x,
                         residual,
                         strong_set,
                         this->x_centers,
                         this->x_scales,
                         Eigen::VectorXd::Ones(n),
                         jit_normalization);

          auto violations = setDiff(
            kktCheck(gradient, beta, lambda_curr, strong_set), working_set);

          if (violations.empty()) {
            updateGradient(gradient,
                           x,
                           residual,
                           full_set,
                           this->x_centers,
                           this->x_scales,
                           Eigen::VectorXd::Ones(n),
                           jit_normalization);

            violations = setDiff(
              kktCheck(gradient, beta, lambda_curr, full_set), working_set);
            if (violations.empty()) {
              break;
            } else {
              working_set = setUnion(working_set, violations);
            }
          } else {
            working_set = setUnion(working_set, violations);
          }
        } else {
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
                  x,
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

    // Store everything for this step of the path
    auto [beta0_out, beta_out] = rescaleCoefficients(
      beta0, beta, this->x_centers, this->x_scales, this->intercept);

    alpha_prev = alpha_curr;

    // Compute early stopping criteria
    double dev = loss->deviance(eta, y);
    double dev_ratio = 1 - dev / null_deviance;
    double dev_change = path_step == 0 ? 1.0 : 1 - dev / dev_prev;
    dev_prev = dev;

    std::vector<std::vector<int>> clusters;

    if (return_clusters) {
      Clusters beta_clusters(beta);
      clusters = beta_clusters.getClusters();
    }

    SlopeFit fit{
      beta0_out, beta_out.sparseView(), clusters,          alpha_curr, lambda,
      dev,       null_deviance,         primals,           duals,      time,
      it,        this->centering_type,  this->scaling_type
    };

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

template<typename T>
SlopeFit
Slope::fit(T& x,
           const Eigen::MatrixXd& y_in,
           const double alpha,
           Eigen::ArrayXd lambda)
{
  Eigen::ArrayXd alpha_arr(1);
  alpha_arr(0) = alpha;
  SlopePath res = path(x, y_in, alpha_arr, lambda);

  return { res(0) };
};

void
Slope::setSolver(const std::string& solver)
{
  validateOption(solver, { "auto", "pgd", "hybrid", "fista" }, "solver");
  this->solver_type = solver;
}

void
Slope::setIntercept(bool intercept)
{
  this->intercept = intercept;
}

void
Slope::setNormalization(const std::string& type)
{
  validateOption(
    type, { "standardization", "min_max", "max_abs", "none" }, "type");

  if (type == "standardization") {
    this->centering_type = "mean";
    this->scaling_type = "sd";
  } else if (type == "none") {
    this->centering_type = "none";
    this->scaling_type = "none";
  } else if (type == "min_max") {
    this->centering_type = "min";
    this->scaling_type = "range";
  } else if (type == "max_abs") {
    this->centering_type = "none";
    this->scaling_type = "max_abs";
  }
}

void
Slope::setCentering(const std::string& type)
{
  validateOption(type, { "mean", "min", "none" }, "type");
  this->centering_type = type;
}

void
Slope::setCentering(const Eigen::VectorXd& x_centers)
{
  this->x_centers = x_centers;
  this->centering_type = "manual";
}

void
Slope::setScaling(const std::string& type)
{
  validateOption(
    type, { "sd", "l1", "l2", "range", "max_abs", "none" }, "type");
  this->scaling_type = type;
}

void
Slope::setScaling(const Eigen::VectorXd& x_scales)
{
  this->x_scales = x_scales;
  this->scaling_type = "manual";
}

void
Slope::setUpdateClusters(bool update_clusters)
{
  this->update_clusters = update_clusters;
}

void
Slope::setReturnClusters(const bool return_clusters)
{
  this->return_clusters = return_clusters;
}

void
Slope::setAlphaType(const std::string& alpha_type)
{
  validateOption(alpha_type, { "path", "estimate" }, "alpha_type");
  this->alpha_type = alpha_type;
}

void
Slope::setAlphaMinRatio(double alpha_min_ratio)
{
  if (alpha_min_ratio <= 0 || alpha_min_ratio >= 1) {
    throw std::invalid_argument("alpha_min_ratio must be in (0, 1)");
  }
  this->alpha_min_ratio = alpha_min_ratio;
}

void
Slope::setLearningRateDecr(double learning_rate_decr)
{
  if (learning_rate_decr <= 0 || learning_rate_decr >= 1) {
    throw std::invalid_argument("learning_rate_decr must be in (0, 1)");
  }
  this->learning_rate_decr = learning_rate_decr;
}

void
Slope::setQ(double q)
{
  if (q < 0 || q > 1) {
    throw std::invalid_argument("q must be between 0 and 1");
  }
  this->q = q;
}

void
Slope::setOscarParameters(const double theta1, const double theta2)
{
  if (theta1 < 0) {
    throw std::invalid_argument("theta1 must be between 0 and 1");
  }

  if (theta2 < 0) {
    throw std::invalid_argument("theta2 must be between 0 and 1");
  }

  this->theta1 = theta1;
  this->theta2 = theta2;
}

void
Slope::setTol(double tol)
{
  if (tol < 0) {
    throw std::invalid_argument("tol must be non-negative");
  }
  this->tol = tol;
}

void
Slope::setMaxIterations(int max_it)
{
  if (max_it < 1) {
    throw std::invalid_argument("max_it_outer must be >= 1");
  }
  this->max_it = max_it;
}

void
Slope::setPathLength(int path_length)
{
  if (path_length < 1) {
    throw std::invalid_argument("path_length must be >= 1");
  }
  this->path_length = path_length;
}

void
Slope::setHybridCdIterations(int cd_iterations)
{
  if (cd_iterations < 0) {
    throw std::invalid_argument("cd_iterations must be >= 0");
  }
  this->cd_iterations = cd_iterations;
}

void
Slope::setLambdaType(const std::string& lambda_type)
{
  validateOption(
    lambda_type, { "bh", "gaussian", "oscar", "lasso" }, "lambda_type");

  this->lambda_type = lambda_type;
}

void
Slope::setLoss(const std::string& loss_type)
{
  validateOption(loss_type,
                 { "quadratic", "logistic", "poisson", "multinomial" },
                 "loss_type");
  this->loss_type = loss_type;
}

void
Slope::setScreening(const std::string& screening_type)
{
  validateOption(screening_type, { "strong", "none" }, "screening_type");
  this->screening_type = screening_type;
}

void
Slope::setModifyX(const bool modify_x)
{
  this->modify_x = modify_x;
}

void
Slope::setDevChangeTol(const double dev_change_tol)
{
  if (dev_change_tol < 0 || dev_change_tol > 1) {
    throw std::invalid_argument("dev_change_tol must be in [0, 1]");
  }

  this->dev_change_tol = dev_change_tol;
}

void
Slope::setDevRatioTol(const double dev_ratio_tol)
{
  if (dev_ratio_tol < 0 || dev_ratio_tol > 1) {
    throw std::invalid_argument("dev_ratio_tol must be in [0, 1]");
  }

  this->dev_ratio_tol = dev_ratio_tol;
}

void
Slope::setMaxClusters(const int max_clusters)
{
  if (max_clusters < 1) {
    throw std::invalid_argument("max_clusters must be >= 1");
  }

  this->max_clusters = max_clusters;
}

void
Slope::setDiagnostics(const bool collect_diagnostics)
{
  this->collect_diagnostics = collect_diagnostics;
}

void
Slope::setAlphaEstimationMaxIterations(const int alpha_est_maxit)
{
  this->alpha_est_maxit = alpha_est_maxit;
}

int
Slope::getAlphaEstimationMaxIterations() const
{
  return alpha_est_maxit;
}

bool
Slope::getFitIntercept() const
{
  return intercept;
}

const std::string&
Slope::getLossType()
{
  return loss_type;
}

/// @cond
// Explicit template instantiations
template SlopePath
Slope::path<Eigen::MatrixXd>(Eigen::MatrixXd&,
                             const Eigen::MatrixXd&,
                             Eigen::ArrayXd,
                             Eigen::ArrayXd);

template SlopePath
Slope::path<Eigen::SparseMatrix<double>>(Eigen::SparseMatrix<double>&,
                                         const Eigen::MatrixXd&,
                                         Eigen::ArrayXd,
                                         Eigen::ArrayXd);

template SlopeFit
Slope::fit<Eigen::MatrixXd>(Eigen::MatrixXd&,
                            const Eigen::MatrixXd&,
                            const double,
                            Eigen::ArrayXd);

template SlopeFit
Slope::fit<Eigen::SparseMatrix<double>>(Eigen::SparseMatrix<double>&,
                                        const Eigen::MatrixXd&,
                                        const double,
                                        Eigen::ArrayXd);
/// @endcond

} // namespace slope
