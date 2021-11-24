#pragma once

#include <RcppArmadillo.h>

struct Results
{
  const arma::cube betas;
  const arma::field<arma::uvec> active_sets;
  const arma::uvec passes;
  const std::vector<std::vector<double>> primals;
  const std::vector<std::vector<double>> duals;
  const std::vector<std::vector<double>> time;
  const arma::uvec n_unique;
  const arma::vec deviance_ratio;
  const double null_deviance;
  const std::vector<std::vector<unsigned>> violations;
  const arma::vec alpha;
  const arma::vec lambda;

  Results(const arma::cube betas,
          const arma::field<arma::uvec> active_sets,
          const arma::uvec passes,
          const std::vector<std::vector<double>> primals,
          const std::vector<std::vector<double>> duals,
          const std::vector<std::vector<double>> time,
          const arma::uvec n_unique,
          const arma::vec deviance_ratio,
          const double null_deviance,
          const std::vector<std::vector<unsigned>> violations,
          const arma::vec alpha,
          const arma::vec lambda);
};
