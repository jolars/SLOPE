#pragma once

#include <RcppArmadillo.h>

arma::uvec
strongSet(const arma::mat& gradient_prev,
          const arma::vec& lambda,
          const arma::vec& lambda_prev,
          const bool intercept);
