#include "Results.h"
#include "SolverResults.h"
#include "families/families.h"
#include "kktCheck.h"
#include "regularizationPath.h"
#include "rescale.h"
#include "screening.h"
#include "standardize.h"
#include <RcppArmadillo.h>
#include <memory>

template<typename T>
Results
SLOPE(T& x,
      arma::mat& y,
      const std::string family_choice,
      const bool intercept,
      arma::vec lambda,
      arma::vec alpha,
      const std::string lambda_type,
      const std::string alpha_type,
      const double alpha_min_ratio,
      const double q,
      const double theta1,
      const double theta2,
      const bool center,
      const std::string scale,
      const arma::rowvec y_center,
      const arma::rowvec y_scale,
      const double tol_dev_ratio,
      const double tol_dev_change,
      const arma::uword max_variables,
      bool screen,
      const std::string screen_alg,
      const std::string solver,
      const arma::uword max_passes,
      const double tol_rel_gap,
      const double tol_infeas,
      const double tol_abs,
      const double tol_rel,
      const double tol_rel_coef_change,
      const int prox_method_choice,
      const bool diagnostics,
      const arma::uword verbosity)
{
  using namespace Rcpp;
  using namespace arma;

  using std::endl;
  using std::setw;
  using std::showpoint;

  // significant digits
  Rcout.precision(4);

  auto n = x.n_rows;
  auto p = x.n_cols;
  auto m = y.n_cols;

  auto prox_method = ProxMethod(prox_method_choice);

  rowvec x_center(p, fill::zeros);
  rowvec x_scale(p, fill::ones);

  standardize(x, x_center, x_scale, intercept, center, scale);

  const uword path_length = alpha.n_elem;
  double alpha_max        = 0;

  regularizationPath(alpha,
                     lambda,
                     alpha_max,
                     x,
                     y,
                     x_scale,
                     y_scale,
                     lambda_type,
                     alpha_type,
                     scale,
                     alpha_min_ratio,
                     q,
                     theta1,
                     theta2,
                     family_choice,
                     intercept);

  auto family = setupFamily(family_choice,
                            intercept,
                            prox_method,
                            diagnostics,
                            max_passes,
                            tol_rel_gap,
                            tol_infeas,
                            tol_abs,
                            tol_rel,
                            tol_rel_coef_change,
                            verbosity);

  cube betas(p, m, path_length, fill::zeros);
  mat beta(p, m, fill::zeros);

  uword n_variables = 0;
  uvec n_unique(path_length);

  mat linear_predictor = x * beta;

  double null_deviance = 2 * family->primal(y, linear_predictor);
  vec deviance_ratios(path_length);
  vec deviances(path_length);
  double deviance_change{ 0 };

  mat beta_prev(p, m, fill::zeros);

  uvec passes(path_length);
  std::vector<std::vector<double>> primals;
  std::vector<std::vector<double>> duals;
  std::vector<std::vector<double>> timings;
  std::vector<unsigned> violations;
  std::vector<std::vector<unsigned>> violation_list;

  mat linear_predictor_prev(n, m);
  mat gradient_prev(p, m);
  mat pseudo_gradient_prev(n, m);

  // sets of active predictors
  field<uvec> active_sets(path_length);
  uvec active_set = regspace<uvec>(0, p - 1);
  uvec strong_set;
  uvec previous_set;
  if (intercept)
    previous_set.insert_rows(0, 1);

  // object for use in ADMM
  double rho = 0.0;
  vec z(p);
  vec u(p);
  vec z_subset(z);
  vec u_subset(u);
  // for gaussian case
  mat xx, L, U;
  vec xTy;
  T x_subset;

  if (family->name() == "gaussian" && solver == "admm") {
    // initialize auxiliary variables
    z.zeros();
    u.zeros();
  }

  bool factorized = false;

  SolverResults solver_res;

  uword k = 0;

  while (k < path_length) {

    violations.clear();

    if (screen) {
      // NOTE(JL): the screening rules should probably not be used if
      // the coefficients from the previous fit are already very dense

      gradient_prev = family->gradient(x, y, x * beta_prev);

      double alpha_prev = (k == 0) ? alpha_max : alpha(k - 1);

      strong_set = strongSet(
        gradient_prev, lambda * alpha(k), lambda * alpha_prev, intercept);

      previous_set = find(any(beta_prev != 0, 1));

      if (intercept)
        previous_set = setUnion(previous_set, { 0 });

      strong_set = setUnion(strong_set, previous_set);

      if (screen_alg == "previous") {
        active_set = previous_set;
      } else {
        active_set = strong_set;
      }
    }

    if (active_set.n_elem == p / m || !screen) {

      // stop screening
      screen = false;

      // all features active
      // factorize once if fitting all
      if (!factorized && family->name() == "gaussian" && solver == "admm") {
        // precompute x^Ty
        xTy = x.t() * y;

        // precompute X^tX or XX^t (if wide) and factorize
        if (n >= p) {
          xx = x.t() * x;
        } else {
          xx = x * x.t();
        }

        vec eigval = eig_sym(xx);

        if (lambda.max() * alpha(k) == 0) {
          rho = eigval.max();
        } else {
          rho = std::pow(eigval.max(), 1 / 3) *
                std::pow(lambda.max() * alpha(k), 2 / 3);
        }

        if (n >= p) {
          xx.diag() += rho;
        } else {
          xx /= rho;
          xx.diag() += 1;
        }

        U = chol(xx);
        L = U.t();

        factorized = true;
      }

      solver_res = family->fit(
        x, y, beta, z, u, L, U, xTy, lambda * alpha(k), rho, solver);
      passes(k) = solver_res.passes;
      beta      = solver_res.beta;

      if (diagnostics) {
        primals.push_back(solver_res.primals);
        duals.push_back(solver_res.duals);
        timings.push_back(solver_res.time);
        violation_list.push_back(violations);
      }

    } else {

      bool kkt_violation = true;

      while (kkt_violation) {

        kkt_violation = false;

        x_subset = matrixSubset(x, active_set);

        if (active_set.n_elem == 0) {
          // null model
          beta.zeros();
          passes(k) = 0;

        } else {

          if (family->name() == "gaussian" && solver == "admm") {
            if (x_subset.n_rows >= x_subset.n_cols) {
              xx = x_subset.t() * x_subset;
            } else {
              xx = x_subset * x_subset.t();
            }

            vec eigval = eig_sym(xx);

            if (lambda.max() * alpha(k) == 0) {
              rho = eigval.max();
            } else {
              rho = std::pow(eigval.max(), 2 / 3) *
                    std::pow(lambda.max() * alpha(k), 1 / 3);
            }

            if (x_subset.n_rows < x_subset.n_cols)
              xx /= rho;

            xx.diag() += rho;

            U = chol(xx);
            L = U.t();

            xTy = x_subset.t() * y;

            z_subset = z(active_set);
            u_subset = u(active_set);
          }

          uword n_active =
            (active_set.n_elem - static_cast<uword>(intercept)) * m;

          solver_res = family->fit(x_subset,
                                   y,
                                   beta.rows(active_set),
                                   z_subset,
                                   u_subset,
                                   L,
                                   U,
                                   xTy,
                                   lambda.head(n_active) * alpha(k),
                                   rho,
                                   solver);

          if (family->name() == "gaussian" && solver == "admm") {
            z(active_set) = z_subset;
            u(active_set) = u_subset;
          }

          beta.rows(active_set) = solver_res.beta;
          passes(k)             = solver_res.passes;
        }

        uvec check_failures;
        uword n_strong =
          (strong_set.n_elem - static_cast<uword>(intercept)) * m;

        if (screen_alg == "previous" && n_strong > 0) {
          // check against strong set
          x_subset = matrixSubset(x, strong_set);
          gradient_prev =
            family->gradient(x_subset, y, x_subset * beta.rows(strong_set));
          uvec tmp = kktCheck(gradient_prev,
                              beta.rows(strong_set),
                              lambda.head(n_strong) * alpha(k),
                              tol_infeas,
                              intercept);

          uvec strong_failures = strong_set(tmp);
          check_failures       = setDiff(strong_failures, active_set);

          kkt_violation = check_failures.n_elem > 0;

          if (diagnostics)
            violations.push_back(check_failures.n_elem);
        }

        if (!kkt_violation) {
          // check against whole set
          gradient_prev = family->gradient(x, y, x * beta);

          uvec tmp = kktCheck(
            gradient_prev, beta, lambda * alpha(k), tol_infeas, intercept);

          check_failures = setDiff(tmp, active_set);

          kkt_violation = check_failures.n_elem > 0;

          if (diagnostics)
            violations.push_back(check_failures.n_elem);
        }

        active_set = setUnion(check_failures, active_set);

        checkUserInterrupt();
      }

      if (diagnostics) {
        primals.push_back(solver_res.primals);
        duals.push_back(solver_res.duals);
        timings.push_back(solver_res.time);
        violation_list.push_back(violations);
      }
    }

    // store coefficients and intercept
    double deviance       = solver_res.deviance;
    double deviance_ratio = 1.0 - deviance / null_deviance;
    deviances(k)          = deviance;
    deviance_ratios(k)    = deviance_ratio;

    deviance_change =
      k == 0 ? 0.0 : std::abs((deviances(k - 1) - deviance) / deviances(k - 1));

    betas.slice(k) = beta;
    beta_prev      = beta;

    active_sets(k) = active_set;
    uword n_coefs  = accu(any(beta != 0, 1));
    n_variables    = n_coefs;
    n_unique(k)    = unique(abs(nonzeros(beta))).eval().n_elem;

    if (verbosity >= 1)
      Rcout << showpoint << "penalty: " << setw(2) << k << ", dev: " << setw(7)
            << deviance << ", dev ratio: " << setw(7) << deviance_ratio
            << ", dev change: " << setw(7) << deviance_change
            << ", n var: " << setw(5) << n_variables
            << ", n unique: " << setw(5) << n_unique(k) << endl;

    if (n_coefs > 0 && k > 0) {
      // stop path if fractional deviance change is small
      if (deviance_change < tol_dev_change || deviance_ratio > tol_dev_ratio) {
        k++;
        break;
      }
    }

    if (n_unique(k) > max_variables) {
      k++;
      break;
    }

    checkUserInterrupt();

    k++;
  }

  betas.resize(p, m, k);
  passes.resize(k);
  alpha.resize(k);
  n_unique.resize(k);
  deviances.resize(k);
  deviance_ratios.resize(k);
  active_sets = active_sets.rows(0, std::max(static_cast<int>(k - 1), 0));

  rescale(betas, x_center, x_scale, y_center, y_scale, intercept);

  // rescale alpha depending on standardization settings
  if (scale == "l2") {
    alpha /= std::sqrt(n);
  } else if (scale == "sd" || scale == "none") {
    alpha /= n;
  }

  return Results{ betas,         active_sets,    passes,   primals,
                  duals,         timings,        n_unique, deviance_ratios,
                  null_deviance, violation_list, alpha,    lambda };
}
