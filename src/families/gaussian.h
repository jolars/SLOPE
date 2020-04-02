#pragma once

#include <RcppArmadillo.h>
#include "family.h"
#include "../results.h"
#include "../utils.h"
#include "../infeasibility.h"
#include "../prox.h"

using namespace Rcpp;
using namespace arma;

class Gaussian : public Family {
private:
  double alpha = 1.5;

public:
  template <typename... Ts>
  Gaussian(Ts... args) : Family(std::forward<Ts>(args)...) {}

  double primal(const mat& y, const mat& lin_pred)
  {
    return 0.5*pow(norm(y - lin_pred), 2);
  }

  double dual(const mat& y, const mat& lin_pred)
  {
    using namespace std;
    return 0.5*pow(norm(y, 2), 2) - 0.5*pow(norm(lin_pred, 2), 2);
  }

  mat pseudoGradient(const mat& y, const mat& lin_pred)
  {
    return lin_pred - y;
  }

  rowvec fitNullModel(const mat& y, const uword n_classes)
  {
    return mean(y);
  }

  std::string name()
  {
    return "gaussian";
  }

  Results ADMM(const mat& x,
               const mat& y,
               mat beta,
               vec& z,
               vec& u,
               const mat& L,
               const mat& U,
               const vec& xTy,
               vec lambda,
               double rho)
  {
    return ADMMImpl(x, y, beta, z, u, L, U, xTy, lambda, rho);
  }

  Results ADMM(const sp_mat& x,
               const mat& y,
               mat beta,
               vec& z,
               vec& u,
               const mat& L,
               const mat& U,
               const vec& xTy,
               vec lambda,
               double rho)
  {
    return ADMMImpl(x, y, beta, z, u, L, U, xTy, lambda, rho);
  }

  // ADMM
  template <typename T>
  Results ADMMImpl(const T& x,
                   const mat& y,
                   mat beta,
                   vec& z,
                   vec& u,
                   const mat& L,
                   const mat& U,
                   const vec& xTy,
                   vec lambda,
                   double rho)
  {
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

      q = xTy + rho*(z - u);

      if (n >= p) {
        beta = solve(trimatu(U), solve(trimatl(L), q));
      } else {
        beta = q/rho -
          (x.t() * solve(trimatu(U), solve(trimatl(L), x*q)))/(rho*rho);
      }

      z_old = z;
      beta_hat = alpha*beta + (1 - alpha)*z_old;

      z = beta_hat + u;
      z.tail(lambda.n_elem) = prox(z.tail(lambda.n_elem), lambda/rho);

      u += (beta_hat - z);

      double r_norm = norm(beta - z);
      double s_norm = norm(rho*(z - z_old));

      double eps_primal =
        std::sqrt(n)*tol_abs + tol_rel*std::max(norm(beta), norm(z));
      double eps_dual = std::sqrt(n)*tol_abs + tol_rel*norm(rho*u);

      if (diagnostics) {
        primals.push_back(r_norm);
        duals.push_back(s_norm);
        time.push_back(timer.toc());
      }

      if (verbosity >= 3) {
        Rcout << "pass: "              << passes
              << ", primal residual: " << r_norm
              << ", dual residual: "   << s_norm
              << std::endl;
      }

      if (r_norm < eps_primal && s_norm < eps_dual)
        break;

      checkUserInterrupt();
    }

    double deviance = 2*primal(y, x*z);

    Results res{z,
                passes,
                primals,
                duals,
                time,
                deviance};

    return res;
  }
};


