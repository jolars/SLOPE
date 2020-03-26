#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

#include "solver.h"
#include "fista.h"

// helper to choose solver
template <typename... Ts>
inline std::unique_ptr<Solver> setupSolver(const std::string& solver_choice,
                                           Ts... args)
{
  // if (family_choice == "fista")
    return std::unique_ptr<FISTA>(new FISTA{std::forward<Ts>(args)...});
  // else if (family_choice == "poisson")
  //   return std::unique_ptr<Poisson>(new Poisson{std::forward<Ts>(args)...});
  // else if (family_choice == "multinomial")
  //   return std::unique_ptr<Multinomial>(new Multinomial{std::forward<Ts>(args)...});
  // else
  //   return std::unique_ptr<Gaussian>(new Gaussian{std::forward<Ts>(args)...});
}
