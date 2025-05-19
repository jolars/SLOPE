#include "slope.h"
#include "solvers/setup_solver.h"
#include "utils.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <set>

namespace slope {

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
Slope::setRelaxTol(double tol)
{
  if (tol < 0) {
    throw std::invalid_argument("tol must be non-negative");
  }

  this->tol_relax = tol;
}

void
Slope::setRelaxMaxOuterIterations(int max_it)
{
  if (max_it < 1) {
    throw std::invalid_argument("max_it must be >= 1");
  }
  this->max_it_outer_relax = max_it;
}

void
Slope::setRelaxMaxInnerIterations(int max_it)
{
  if (max_it < 1) {
    throw std::invalid_argument("max_it_outer must be >= 1");
  }
  this->max_it_inner_relax = max_it;
}

void
Slope::setMaxIterations(int max_it)
{
  if (max_it < 1) {
    throw std::invalid_argument("max_it must be >= 1");
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

} // namespace slope
