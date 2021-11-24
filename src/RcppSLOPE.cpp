#include "SLOPE.h"
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List
sparseSLOPE(arma::sp_mat x, arma::mat y, const Rcpp::List control)
{
  return SLOPE(x, y, control);
}

// [[Rcpp::export]]
Rcpp::List
denseSLOPE(arma::mat x, arma::mat y, const Rcpp::List control)
{
  return SLOPE(x, y, control);
}
