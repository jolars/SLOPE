#pragma once

#include <RcppArmadillo.h>
#include "../results.h"
#include "../utils.h"
#include "../infeasibility.h"
#include "../prox.h"

using namespace Rcpp;
using namespace arma;

class Solver {
protected:
  const bool intercept;
  const bool diagnostics;
  const uword max_passes;
  const double tol_rel_gap;
  const double tol_infeas;
  const double tol_abs;
  const double tol_rel;
  const uword verbosity;

public:
  Solver(const bool intercept,
         const bool diagnostics,
         const uword max_passes,
         const double tol_rel_gap,
         const double tol_infeas,
         const double tol_abs,
         const double tol_rel,
         const uword verbosity)
    : intercept(intercept),
      diagnostics(diagnostics),
      max_passes(max_passes),
      tol_rel_gap(tol_rel_gap),
      tol_infeas(tol_infeas),
      tol_abs(tol_abs),
      tol_rel(tol_rel),
      verbosity(verbosity) {}

  virtual Results fit(const mat& x,
                      const mat& y,
                      mat beta,
                      vec& z,
                      vec& u,
                      const mat& L,
                      const mat& U,
                      const vec& xTy,
                      vec lambda,
                      double rho) = 0;
};
