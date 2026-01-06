#include "setup_model.h"
#include "slope/clusters.h"
#include "slope/slope.h"
#include "slope/threads.h"
#include <RcppEigen.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

static void
chkIntFn(void* dummy)
{
  R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of
// your context
bool
checkInterrupt()
{
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

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

  slope::SlopePath path = model.path(x, y, alpha, lambda, checkInterrupt);

  std::vector<Eigen::SparseMatrix<double>> patterns;

  bool return_patterns = Rcpp::as<bool>(control["patterns"]);

  auto coefs = path.getCoefs(false);

  if (return_patterns) {
    for (auto coef : coefs) {
      Eigen::SparseMatrix<double> U = slope::patternMatrix(coef).cast<double>();
      patterns.push_back(U);
    }
  }

  Eigen::SparseMatrix<double> pattern_matrix;

  if (gamma < 1) {
    path = model.relax(path, x, y, gamma);
  }

  return Rcpp::List::create(
    Named("intercepts") = wrap(path.getIntercepts()),
    Named("intercepts_scaled") = wrap(path.getIntercepts(false)),
    Named("betas") = wrap(path.getCoefs()),
    Named("betas_scaled") = wrap(coefs),
    Named("passes") = wrap(path.getPasses()),
    Named("primals") = wrap(path.getPrimals()),
    Named("duals") = wrap(path.getDuals()),
    Named("time") = wrap(path.getTime()),
    Named("deviance_ratio") = wrap(path.getDevianceRatios()),
    Named("null_deviance") = wrap(path.getNullDeviance()),
    Named("alpha") = wrap(path.getAlpha()),
    Named("lambda") = wrap(path.getLambda()),
    Named("patterns") = wrap(patterns));
}

// [[Rcpp::export]]
Rcpp::List
sparseSLOPE(Eigen::SparseMatrix<double>& x,
            const Eigen::MatrixXd& y,
            const Rcpp::List& control)
{
  if (!y.allFinite()) {
    throw std::invalid_argument("Response must not contain missing values");
  }
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
  if (!y.allFinite()) {
    throw std::invalid_argument("Response must not contain missing values");
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

  if (!y.allFinite()) {
    throw std::invalid_argument("Response must not contain missing values");
  }

  return callSLOPE(x_map, y, control);
}
