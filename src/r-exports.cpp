#include <RcppArmadillo.h>
#include "prox.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat sorted_l1_prox(const arma::mat& x, const arma::vec& lambda)
{
  return prox(x, lambda);
}


