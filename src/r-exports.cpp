#include "slope/regularization_sequence.h"
#include "slope/sorted_l1_norm.h"
#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd
sortedL1ProxCpp(const Eigen::MatrixXd& x, const Eigen::ArrayXd& lambda)
{
  slope::SortedL1Norm norm;

  return norm.prox(x, lambda);
}

// [[Rcpp::export]]
Eigen::ArrayXd
lambdaSequenceCpp(const int n_lambda,
                  const double q,
                  const double theta1,
                  const double theta2,
                  const std::string& lambda_type,
                  const int n)
{
  return slope::lambdaSequence(n_lambda, q, lambda_type, n, theta1, theta2);
}
