#pragma once

#include <RcppArmadillo.h>

struct Results
{
  arma::mat beta;
  arma::uword passes;
  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> time;
  double deviance;

  Results();

  Results(arma::mat beta,
          arma::uword passes,
          std::vector<double> primals,
          std::vector<double> duals,
          std::vector<double> time,
          double deviance);
};
