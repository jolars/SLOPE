#pragma once

#include <RcppArmadillo.h>

void
standardize(arma::mat& x,
            arma::rowvec& x_center,
            arma::rowvec& x_scale,
            bool intercept,
            bool center,
            std::string scale);

void
standardize(arma::sp_mat& x,
            arma::rowvec& x_center,
            arma::rowvec& x_scale,
            bool intercept,
            bool center,
            std::string scale);
