#include "lambdaSequence.h"
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::vec
lambdaSequence(const arma::uword n_lambda,
               const double q,
               const double theta1,
               const double theta2,
               const std::string lambda_type,
               const arma::uword n)
{
  using namespace arma;

  vec lambda(n_lambda);

  if (lambda_type == "gaussian" || lambda_type == "bh") {
    lambda = regspace(1, n_lambda) * q / (2 * n_lambda);

    lambda.transform(
      [](double val) { return Rf_qnorm5(1.0 - val, 0.0, 1.0, 1, 0); });

    if (lambda_type == "gaussian" && n_lambda > 1) {
      double sum_sq = 0.0;

      for (uword i = 1; i < n_lambda; ++i) {
        sum_sq += std::pow(lambda(i - 1), 2);
        double w = std::max(1.0, static_cast<double>(n) - i - 1.0);
        lambda(i) *= std::sqrt(1.0 + sum_sq / w);
      }

      // ensure non-increasing lambda
      for (arma::uword i = 1; i < n_lambda; ++i) {
        if (lambda(i - 1) < lambda(i)) {
          lambda(i) = lambda(i - 1);
        }
      }
    }
  } else if (lambda_type == "oscar") {
    lambda = theta1 + theta2 * (n_lambda - regspace(1, n_lambda));
  } else if (lambda_type == "lasso") {
    lambda.ones();
  }

  return lambda;
}
