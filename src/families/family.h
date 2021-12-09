#pragma once

#include "../infeasibility.h"
#include "../prox.h"
#include "../SolverResults.h"
#include "../utils.h"
#include <RcppArmadillo.h>

class Family
{
protected:
  const bool intercept;
  const ProxMethod prox_method;
  const bool diagnostics;
  const arma::uword max_passes;
  const double tol_rel_gap;
  const double tol_infeas;
  const double tol_abs;
  const double tol_rel;
  const double tol_rel_coef_change;
  const arma::uword verbosity;

public:
  Family(const bool intercept,
         const ProxMethod prox_method,
         const bool diagnostics,
         const arma::uword max_passes,
         const double tol_rel_gap,
         const double tol_infeas,
         const double tol_abs,
         const double tol_rel,
         const double tol_rel_coef_change,
         const arma::uword verbosity)
    : intercept(intercept)
    , prox_method(prox_method)
    , diagnostics(diagnostics)
    , max_passes(max_passes)
    , tol_rel_gap(tol_rel_gap)
    , tol_infeas(tol_infeas)
    , tol_abs(tol_abs)
    , tol_rel(tol_rel)
    , tol_rel_coef_change(tol_rel_coef_change)
    , verbosity(verbosity)
  {}

  virtual ~Family(){}

  virtual double primal(const arma::mat& y, const arma::mat& lin_pred) = 0;

  virtual double dual(const arma::mat& y, const arma::mat& lin_pred) = 0;

  // this is not really the true gradient; it needs to multiplied by X^T
  virtual arma::mat pseudoGradient(const arma::mat& y,
                                   const arma::mat& lin_pred) = 0;

  template<typename T>
  arma::mat gradient(const T& x, const arma::mat& y, const arma::mat& lin_pred)
  {
    return x.t() * pseudoGradient(y, lin_pred);
  }

  virtual arma::rowvec fitNullModel(const arma::mat& y,
                                    const arma::uword n_classes) = 0;

  virtual std::string name() = 0;

  template<typename T>
  SolverResults fit(const T& x,
                    const arma::mat& y,
                    arma::mat beta,
                    arma::vec& z,
                    arma::vec& u,
                    const arma::mat& L,
                    const arma::mat& U,
                    const arma::vec& xTy,
                    arma::vec lambda,
                    double rho,
                    const std::string solver)
  {
    if (solver == "admm")
      return ADMM(x, y, beta, z, u, L, U, xTy, lambda, rho);
    else
      return FISTA(x, y, beta, lambda);
  }

  // FISTA implementation
  template<typename T>
  SolverResults FISTAImpl(const T& x,
                          const arma::mat& y,
                          arma::mat beta,
                          arma::vec lambda)
  {
    using namespace arma;
    using namespace Rcpp;

    uword n      = y.n_rows;
    uword p      = x.n_cols;
    uword m      = beta.n_cols;
    uword pmi    = lambda.n_elem;
    uword p_rows = pmi / m;

    mat beta_tilde(beta);
    mat beta_tilde_old(beta);
    mat beta_prev(beta);

    mat lin_pred(n, m);
    mat grad(p, m, fill::zeros);

    double learning_rate = 1.0;

    // line search parameters
    double eta = 0.5;

    // FISTA parameters
    double t = 1;

    // diagnostics
    wall_clock timer;
    std::vector<double> primals;
    std::vector<double> duals;
    std::vector<double> time;

    if (diagnostics) {
      primals.reserve(max_passes);
      duals.reserve(max_passes);
      time.reserve(max_passes);
      timer.tic();
    }

    // main loop
    uword passes = 0;
    while (passes < max_passes) {
      lin_pred = x * beta;

      double g = primal(y, lin_pred);

      mat tmp = beta;

      if (intercept)
        tmp.shed_row(0);

      vec tmp_vectorized = sort(abs(vectorise(tmp)), "descending");

      double h = dot(tmp_vectorized, lambda);
      double f = g + h;
      double G = dual(y, lin_pred);

      grad = gradient(x, y, lin_pred);
      tmp  = grad;

      if (intercept)
        tmp.shed_row(0);

      tmp_vectorized = vectorise(tmp);

      double infeas =
        lambda.n_elem > 0.0 ? infeasibility(tmp_vectorized, lambda) : 0.0;

      if (verbosity >= 3) {
        Rcout << "pass: " << passes
              << ", duality-gap: " << std::abs(f - G) / std::abs(f)
              << ", infeasibility: " << infeas << std::endl;
      }

      double small = std::sqrt(datum::eps);

      bool optimal =
        (std::abs(f - G) / std::max(small, std::abs(f)) < tol_rel_gap);

      bool feasible = lambda.n_elem > 0.0
                        ? infeas <= std::max(small, tol_infeas * lambda(0))
                        : true;

      // check change in coefficients
      double max_change = abs(vectorise(beta - beta_prev)).max();
      double max_size   = abs(vectorise(beta)).max();

      bool all_zero = (max_size == 0.0) && (max_change == 0.0);
      bool small_change =
        (max_size != 0.0) && (max_change / max_size <= tol_rel_coef_change);

      bool small_coef_change = all_zero || small_change;

      if (diagnostics) {
        time.push_back(timer.toc());
        primals.push_back(f);
        duals.push_back(G);
      }

      if (optimal && feasible && passes > 0 && small_coef_change)
        break;

      beta_tilde_old = beta_tilde;

      double g_old = g;
      double t_old = t;

      // Backtracking line search
      while (true) {
        // Update coefficients
        beta_tilde = beta - learning_rate * grad;

        if (intercept) {
          mat tmp = beta_tilde;
          tmp.shed_row(0);
          beta_tilde.tail_rows(p_rows) =
            prox(tmp, lambda * learning_rate, prox_method);
        } else {
          beta_tilde = prox(beta_tilde, lambda * learning_rate, prox_method);
        }

        vec d = vectorise(beta_tilde - beta);

        lin_pred = x * beta_tilde;

        g = primal(y, lin_pred);

        double q = g_old + dot(d, vectorise(grad)) +
                   (1.0 / (2 * learning_rate)) * accu(square(d));

        if (q >= g * (1 - 1e-12)) {
          break;
        } else {
          learning_rate *= eta;
        }

        checkUserInterrupt();
      }

      // FISTA step
      t    = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t_old * t_old));
      beta = beta_tilde + (t_old - 1.0) / t * (beta_tilde - beta_tilde_old);

      beta_prev = beta;

      if (passes % 100 == 0)
        checkUserInterrupt();

      ++passes;
    }

    double deviance = 2 * primal(y, lin_pred);

    SolverResults res{ beta, passes, primals, duals, time, deviance };

    return res;
  }

  virtual SolverResults FISTA(const arma::sp_mat& x,
                              const arma::mat& y,
                              arma::mat beta,
                              arma::vec lambda)
  {
    return FISTAImpl(x, y, beta, lambda);
  }

  virtual SolverResults FISTA(const arma::mat& x,
                              const arma::mat& y,
                              arma::mat beta,
                              arma::vec lambda)
  {
    return FISTAImpl(x, y, beta, lambda);
  }

  virtual SolverResults ADMM(const arma::sp_mat& x,
                             const arma::mat& y,
                             arma::mat beta,
                             arma::vec& z,
                             arma::vec& u,
                             const arma::mat& L,
                             const arma::mat& U,
                             const arma::vec& xTy,
                             arma::vec lambda,
                             double rho)
  {
    return ADMMImpl(x, y, beta, z, u, L, U, xTy, lambda, rho);
  }

  virtual SolverResults ADMM(const arma::mat& x,
                             const arma::mat& y,
                             arma::mat beta,
                             arma::vec& z,
                             arma::vec& u,
                             const arma::mat& L,
                             const arma::mat& U,
                             const arma::vec& xTy,
                             arma::vec lambda,
                             double rho)
  {
    return ADMMImpl(x, y, beta, z, u, L, U, xTy, lambda, rho);
  }

  // ADMM implementation
  template<typename T>
  SolverResults ADMMImpl(const T& x,
                         const arma::mat& y,
                         arma::mat beta,
                         arma::vec& z,
                         arma::vec& u,
                         const arma::mat& L,
                         const arma::mat& U,
                         const arma::vec& xTy,
                         arma::vec lambda,
                         double rho)
  {
    Rcpp::stop("ADMM solver is not implemented for this family");

    SolverResults res{};

    return res;
  }
};
