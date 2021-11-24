#pragma once

#include <RcppArmadillo.h>

enum class ProxMethod
{
  stack = 0,
  pava  = 1
};

arma::mat
prox(const arma::mat& beta,
     const arma::vec& lambda,
     const ProxMethod prox_method);
