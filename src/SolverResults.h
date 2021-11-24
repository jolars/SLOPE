#pragma once

#include <RcppArmadillo.h>

struct SolverResults
{
  arma::mat beta;
  arma::uword passes;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> time;
  double deviance;

  SolverResults();

  SolverResults(arma::mat beta,
                arma::uword passes,
                std::vector<double> primals,
                std::vector<double> duals,
                std::vector<double> time,
                double deviance);
};
