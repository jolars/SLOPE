#include "Results.h"
#include "SLOPE.h"
#include <RcppArmadillo.h>

template<typename T>
Rcpp::List
callSLOPE(T& x, arma::mat& y, const Rcpp::List control)
{
  using namespace arma;

  using Rcpp::as;
  using Rcpp::Named;
  using Rcpp::wrap;

  auto family_choice = as<std::string>(control["family"]);
  auto intercept = as<bool>(control["fit_intercept"]);
  auto lambda = as<vec>(control["lambda"]);
  auto alpha = as<vec>(control["alpha"]);
  auto lambda_type = as<std::string>(control["lambda_type"]);
  auto alpha_type = as<std::string>(control["alpha_type"]);
  auto alpha_min_ratio = as<double>(control["alpha_min_ratio"]);
  auto q = as<double>(control["q"]);
  auto theta1 = as<double>(control["theta1"]);
  auto theta2 = as<double>(control["theta2"]);
  auto center = as<bool>(control["center"]);
  auto scale = as<std::string>(control["scale"]);
  auto y_center = as<rowvec>(control["y_center"]);
  auto y_scale = as<rowvec>(control["y_scale"]);
  auto tol_dev_ratio = as<double>(control["tol_dev_ratio"]);
  auto tol_dev_change = as<double>(control["tol_dev_change"]);
  auto max_variables = as<uword>(control["max_variables"]);
  auto screen = as<bool>(control["screen"]);
  auto screen_alg = as<std::string>(control["screen_alg"]);
  auto solver = as<std::string>(control["solver"]);
  auto max_passes = as<uword>(control["max_passes"]);
  auto tol_rel_gap = as<double>(control["tol_rel_gap"]);
  auto tol_infeas = as<double>(control["tol_infeas"]);
  auto tol_abs = as<double>(control["tol_abs"]);
  auto tol_rel = as<double>(control["tol_rel"]);
  auto tol_rel_coef_change = as<double>(control["tol_rel_coef_change"]);
  auto prox_method_choice = as<int>(control["prox_method_choice"]);
  auto diagnostics = as<bool>(control["diagnostics"]);
  auto verbosity = as<uword>(control["verbosity"]);

  Results res = SLOPE(x,
                      y,
                      family_choice,
                      intercept,
                      lambda,
                      alpha,
                      lambda_type,
                      alpha_type,
                      alpha_min_ratio,
                      q,
                      theta1,
                      theta2,
                      center,
                      scale,
                      y_center,
                      y_scale,
                      tol_dev_ratio,
                      tol_dev_change,
                      max_variables,
                      screen,
                      screen_alg,
                      solver,
                      max_passes,
                      tol_rel_gap,
                      tol_infeas,
                      tol_abs,
                      tol_rel,
                      tol_rel_coef_change,
                      prox_method_choice,
                      diagnostics,
                      verbosity);

  return Rcpp::List::create(Named("betas") = wrap(res.betas),
                            Named("active_sets") = wrap(res.active_sets),
                            Named("passes") = wrap(res.passes),
                            Named("primals") = wrap(res.primals),
                            Named("duals") = wrap(res.duals),
                            Named("time") = wrap(res.time),
                            Named("n_unique") = wrap(res.n_unique),
                            Named("violations") = wrap(res.violations),
                            Named("deviance_ratio") = wrap(res.deviance_ratio),
                            Named("null_deviance") = wrap(res.null_deviance),
                            Named("alpha") = wrap(res.alpha),
                            Named("lambda") = wrap(res.lambda));
}

// [[Rcpp::export]]
Rcpp::List
sparseSLOPE(arma::sp_mat x, arma::mat y, const Rcpp::List control)
{
  return callSLOPE(x, y, control);
}

// [[Rcpp::export]]
Rcpp::List
denseSLOPE(arma::mat x, arma::mat y, const Rcpp::List control)
{
  return callSLOPE(x, y, control);
}
