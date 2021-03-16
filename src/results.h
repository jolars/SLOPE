#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

struct Results
{
  mat beta;
  uword passes;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> time;
  double deviance;

  Results() {}

  Results(mat beta,
          uword passes,
          std::vector<double> primals,
          std::vector<double> duals,
          std::vector<double> time,
          double deviance)
    : beta(beta)
    , passes(passes)
    , primals(primals)
    , duals(duals)
    , time(time)
    , deviance(deviance)
  {}
};
