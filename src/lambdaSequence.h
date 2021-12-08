#pragma once

#include <RcppArmadillo.h>

arma::vec
lambdaSequence(const arma::uword n_lambda,
               const double q,
               const double theta1,
               const double theta2,
               const std::string lambda_type,
               const arma::uword n);
