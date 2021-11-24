#pragma once

#include "lambdaMax.h"
#include "lambdaSequence.h"
#include <RcppArmadillo.h>

template<typename T>
void
regularizationPath(arma::vec& alpha,
                   arma::vec& lambda,
                   double& alpha_max,
                   const T& x,
                   const arma::mat& y,
                   const arma::rowvec& x_scale,
                   const arma::rowvec& y_scale,
                   const std::string lambda_type,
                   const std::string alpha_type,
                   const std::string scale,
                   const double alpha_min_ratio,
                   const double q,
                   const double theta1,
                   const double theta2,
                   const std::string family,
                   const bool intercept)
{
  using namespace arma;
  using namespace Rcpp;

  const sword n = x.n_rows;
  const uword m = y.n_cols;

  const sword n_lambda    = lambda.n_elem;
  const uword path_length = alpha.n_elem;

  if (lambda_type != "user") {
    lambda = lambdaSequence(n_lambda, q, theta1, theta2, lambda_type, n);
  }

  vec lambda_max = lambdaMax(x, y, y_scale, m, family, intercept);

  alpha_max =
    (cumsum(sort(abs(lambda_max), "descending")) / cumsum(lambda)).max();

  if (alpha_type == "auto") {
    alpha = exp(
      linspace(log(alpha_max), log(alpha_max * alpha_min_ratio), path_length));
  } else if (alpha_type == "user") {
    // scale alpha to make penalty invariant to number of observations
    if (scale == "l2") {
      alpha *= std::sqrt(n);
    } else if (scale == "sd" || scale == "none") {
      alpha *= n;
    }
  }
}
