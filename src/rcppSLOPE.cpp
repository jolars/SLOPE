#include "slope/slope.h"
#include <RcppEigen.h>

template<typename T>
Rcpp::List
callSLOPE(T& x, const Eigen::MatrixXd& y, const Rcpp::List& control)
{
  using Eigen::ArrayXd;
  using Eigen::VectorXd;
  using Rcpp::as;
  using Rcpp::Named;
  using Rcpp::wrap;

  auto diagnostics = as<bool>(control["diagnostics"]);
  auto alpha = as<ArrayXd>(control["alpha"]);
  auto alpha_min_ratio = as<double>(control["alpha_min_ratio"]);
  auto centering_type = as<std::string>(control["centering_type"]);
  auto centers = as<VectorXd>(control["centers"]);
  auto intercept = as<bool>(control["fit_intercept"]);
  auto lambda = as<ArrayXd>(control["lambda"]);
  auto lambda_type = as<std::string>(control["lambda_type"]);
  auto max_clusters = as<int>(control["max_variables"]);
  auto max_it = as<int>(control["max_passes"]);
  auto loss_type = as<std::string>(control["family"]);
  auto path_length = as<int>(control["path_length"]);
  auto q = as<double>(control["q"]);
  auto scales = as<VectorXd>(control["scales"]);
  auto scaling_type = as<std::string>(control["scaling_type"]);
  auto solver_type = as<std::string>(control["solver"]);
  auto theta1 = as<double>(control["theta1"]);
  auto theta2 = as<double>(control["theta2"]);
  auto tol = as<double>(control["tol"]);
  auto tol_dev_change = as<double>(control["tol_dev_change"]);
  auto tol_dev_ratio = as<double>(control["tol_dev_ratio"]);

  slope::Slope model;

  // Map family to loss_type
  if (loss_type == "gaussian") {
    loss_type = "quadratic";
  } else if (loss_type == "binomial") {
    loss_type = "logistic";
  }

  // overwrite solver type for Poisson regression,
  // since the hybrid solver has convergence issues for this
  // objective
  if (solver_type == "auto") {
    if (loss_type == "poisson") {
      solver_type = "fista";
    }
  }

  model.setAlphaMinRatio(alpha_min_ratio);
  model.setDevChangeTol(tol_dev_change);
  model.setDevRatioTol(tol_dev_ratio);
  model.setIntercept(intercept);
  model.setMaxClusters(max_clusters);
  model.setMaxIterations(max_it);
  model.setLoss(loss_type);
  model.setOscarParameters(theta1, theta2);
  model.setPathLength(path_length);
  model.setQ(q);
  model.setSolver(solver_type);
  model.setTol(tol);
  model.setDiagnostics(diagnostics);

  if (centering_type == "manual") {
    model.setCentering(centers);
  } else {
    model.setCentering(centering_type);
  }

  if (lambda_type != "user") {
    model.setLambdaType(lambda_type);
  }

  if (scaling_type == "manual") {
    model.setScaling(scales);
  } else {
    model.setScaling(scaling_type);
  }

  slope::SlopePath fit = model.path(x, y, alpha, lambda);

  return Rcpp::List::create(
    Named("intercepts") = wrap(fit.getIntercepts()),
    Named("betas") = wrap(fit.getCoefs()),
    Named("passes") = wrap(fit.getPasses()),
    Named("primals") = wrap(fit.getPrimals()),
    Named("duals") = wrap(fit.getDuals()),
    Named("time") = wrap(fit.getTime()),
    Named("deviance_ratio") = wrap(fit.getDevianceRatios()),
    Named("null_deviance") = wrap(fit.getNullDeviance()),
    Named("alpha") = wrap(fit.getAlpha()),
    Named("lambda") = wrap(fit.getLambda()));
}

// [[Rcpp::export]]
Rcpp::List
sparseSLOPE(Eigen::SparseMatrix<double>& x,
            const Eigen::MatrixXd& y,
            const Rcpp::List& control)
{
  return callSLOPE(x, y, control);
}

// [[Rcpp::export]]
Rcpp::List
denseSLOPE(Eigen::MatrixXd& x,
           const Eigen::MatrixXd& y,
           const Rcpp::List& control)
{
  if (!x.allFinite()) {
    throw std::invalid_argument("Input matrix must not contain missing values");
  }
  return callSLOPE(x, y, control);
}
