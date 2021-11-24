#pragma once

#include <RcppArmadillo.h>

void
rescale(arma::cube& betas,
        const arma::rowvec& x_center,
        const arma::rowvec& x_scale,
        const arma::rowvec& y_center,
        const arma::rowvec& y_scale,
        const bool intercept);

