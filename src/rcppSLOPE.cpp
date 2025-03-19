#include "setup_model.h"
#include "slope/slope.h"
#include "slope/threads.h"
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

  slope::Slope model = setupModel(control);

  auto alpha = as<ArrayXd>(control["alpha"]);
  auto lambda = as<ArrayXd>(control["lambda"]);
  auto lambda_type = as<std::string>(control["lambda_type"]);
  auto threads = as<int>(control["threads"]);

  if (lambda_type != "user") {
    model.setLambdaType(lambda_type);
  }

  if (threads > 0) {
    slope::Threads::set(threads);
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
