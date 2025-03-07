#pragma once

#include "slope/slope.h"
#include <RcppEigen.h>

slope::Slope
setupModel(const Rcpp::List& control);
