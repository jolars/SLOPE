#pragma once

#include "../SolverResults.h"
#include "../infeasibility.h"
#include "../prox.h"
#include "../utils.h"
#include "family.h"
#include <RcppArmadillo.h>

class Gaussian : public Family
{
private:
  double alpha = 1.5;

public:
  template<typename... Ts>
  Gaussian(Ts... args)
    : Family(std::forward<Ts>(args)...)
  {}

  double primal(const arma::mat& y, const arma::mat& lin_pred)
  {
    return 0.5 * std::pow(arma::norm(y - lin_pred), 2);
  }

  double dual(const arma::mat& y, const arma::mat& lin_pred)
  {
    return 0.5 * std::pow(arma::norm(y, 2), 2) -
           0.5 * std::pow(arma::norm(lin_pred, 2), 2);
  }

  arma::mat partialGradient(const arma::mat& y, const arma::mat& lin_pred)
  {
    return lin_pred - y;
  }

  arma::rowvec fitNullModel(const arma::mat& y, const arma::uword n_classes)
  {
    return arma::mean(y);
  }

  std::string name() { return "gaussian"; }

  SolverResults ADMM(const arma::mat& x,
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

  SolverResults ADMM(const arma::sp_mat& x,
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

  // ADMM
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
    using namespace arma;
    using namespace Rcpp;

    std::vector<double> primals;
    std::vector<double> duals;
    std::vector<double> infeasibilities;
    std::vector<double> time;

    uword p = x.n_cols;
    uword n = x.n_rows;

    wall_clock timer;

    if (diagnostics)
      timer.tic();

    vec beta_hat(beta);

    vec z_old(z);
    vec q;

    vec lin_pred(n);
    vec grad(p);

    // ADMM loop
    uword passes = 0;

    while (passes < max_passes) {
      ++passes;

      q = xTy + rho * (z - u);

      if (n >= p) {
        beta = solve(trimatu(U), solve(trimatl(L), q));
      } else {
        beta = q / rho - (x.t() * solve(trimatu(U), solve(trimatl(L), x * q))) /
                           (rho * rho);
      }

      z_old    = z;
      beta_hat = alpha * beta + (1 - alpha) * z_old;

      z = beta_hat + u;
      z.tail(lambda.n_elem) =
        prox(z.tail(lambda.n_elem), lambda / rho, prox_method);

      u += (beta_hat - z);

      double r_norm = norm(beta - z);
      double s_norm = norm(rho * (z - z_old));

      double eps_primal =
        std::sqrt(p) * tol_abs + tol_rel * std::max(norm(beta), norm(z));
      double eps_dual = std::sqrt(p) * tol_abs + tol_rel * norm(rho * u);

      if (diagnostics) {
        primals.push_back(r_norm);
        duals.push_back(s_norm);
        time.push_back(timer.toc());
      }

      if (verbosity >= 3) {
        Rcout << "pass: " << passes << ", primal residual: " << r_norm
              << ", dual residual: " << s_norm << std::endl;
      }

      if (r_norm < eps_primal && s_norm < eps_dual)
        break;

      checkUserInterrupt();
    }

    double deviance = 2 * primal(y, x * z);

    SolverResults res{ z, passes, primals, duals, time, deviance };

    return res;
  }
};
