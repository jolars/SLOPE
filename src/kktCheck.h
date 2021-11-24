#pragma once

#include <RcppArmadillo.h>

arma::uvec
kktCheck(arma::mat gradient,
         arma::mat beta,
         const arma::vec& lambda,
         const double tol,
         const bool intercept);
