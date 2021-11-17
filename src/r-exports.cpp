#include "prox.h"
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat
sorted_l1_prox(const arma::mat& x, const arma::vec& lambda)
{
  return prox(x, lambda);
}
