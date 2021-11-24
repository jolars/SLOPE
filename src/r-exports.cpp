#include "prox.h"
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat
sortedL1ProxCpp(const arma::mat& x, const arma::vec& lambda, const int method)
{
  auto prox_method = ProxMethod(method);

  return prox(x, lambda, prox_method);
}
