#include "setup_model.h"
#include "slope/slope.h"
#include "slope/threads.h"
#include <RcppEigen.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

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
  auto gamma = as<double>(control["gamma"]);

  if (lambda_type != "user") {
    model.setLambdaType(lambda_type);
  }

  if (threads > 0) {
    slope::Threads::set(threads);
  }

  slope::SlopePath fit = model.path(x, y, alpha, lambda);

  std::vector<Eigen::SparseMatrix<double>> patterns;

  bool return_patterns = Rcpp::as<bool>(control["patterns"]);

  if (return_patterns) {
    auto clusters = fit.getClusters();

    for (auto cluster : clusters) {
      patterns.push_back(cluster.patternMatrix().cast<double>());
    }
  }

  Eigen::SparseMatrix<double> pattern_matrix;

  if (gamma < 1) {
    fit = model.relax(fit, x, y, gamma);
  }

  return Rcpp::List::create(
    Named("intercepts") = wrap(fit.getIntercepts()),
    Named("intercepts_scaled") = wrap(fit.getIntercepts(false)),
    Named("betas") = wrap(fit.getCoefs()),
    Named("betas_scaled") = wrap(fit.getCoefs(false)),
    Named("passes") = wrap(fit.getPasses()),
    Named("primals") = wrap(fit.getPrimals()),
    Named("duals") = wrap(fit.getDuals()),
    Named("time") = wrap(fit.getTime()),
    Named("deviance_ratio") = wrap(fit.getDevianceRatios()),
    Named("null_deviance") = wrap(fit.getNullDeviance()),
    Named("alpha") = wrap(fit.getAlpha()),
    Named("lambda") = wrap(fit.getLambda()),
    Named("patterns") = wrap(patterns));
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

// [[Rcpp::export]]
Rcpp::List
bigSLOPE(SEXP x, const Eigen::MatrixXd& y, const Rcpp::List& control)
{
  using Eigen::Map;
  using Eigen::MatrixXd;

  Rcpp::XPtr<BigMatrix> bm_ptr(x);

  Map<MatrixXd> x_map =
    Map<MatrixXd>((double*)bm_ptr->matrix(), bm_ptr->nrow(), bm_ptr->ncol());

  return callSLOPE(x_map, y, control);
}
