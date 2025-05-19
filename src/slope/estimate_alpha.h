/**
 * @file
 * @brief Functions for estimating noise level and regularization parameter
 * alpha
 */

#pragma once

#include "ols.h"

namespace slope {

/**
 * @brief Estimates noise (standard error) in a linear model using OLS residuals
 *
 * Calculates the standard error of the residuals from an Ordinary Least Squares
 * (OLS) regression. This is computed as sqrt(sum(residuals²) / (n - p +
 * intercept_term)), where n is the number of observations and p is the number
 * of predictors.
 *
 * @tparam MatrixType The type of matrix used to store the design matrix
 * @param x Design matrix with n observations and p predictors
 * @param y Response matrix
 * @param fit_intercept Whether to include an intercept term in the model
 * @return The estimated noise level (standard error of residuals)
 */
template<typename T>
double
estimateNoise(Eigen::EigenBase<T>& x,
              Eigen::MatrixXd& y,
              const bool fit_intercept)
{
  int n = x.rows();
  int p = x.cols();

  Eigen::VectorXd residuals;
  int df; // Degrees of freedom

  if (p == 0) {
    // Handle the case when X has zero columns
    residuals = y;
    df = n - 1;

    if (fit_intercept) {
      y.array() -= y.mean();
      df -= 1;
    }
  } else {
    // Normal case with predictors
    auto [ols_intercept, ols_coefs] = fitOls(x.derived(), y, fit_intercept);
    residuals = y - x.derived() * ols_coefs;

    if (fit_intercept) {
      residuals.array() -= ols_intercept;
    }

    df = n - p - static_cast<int>(fit_intercept) - 1;
  }

  assert(df > 0);

  return residuals.norm() / std::sqrt(df);
}

} // namespace slope
