/**
 * @file
 * @brief Functions for estimating noise level and regularization parameter
 * alpha
 */

#pragma once

#include "slope.h"
#include "slope/logger.h"
#include "slope/ols.h"
#include "utils.h"

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
template<typename MatrixType>
double
estimateNoise(MatrixType& x, Eigen::MatrixXd& y, const bool fit_intercept)
{
  int n = x.rows();
  int p = x.cols();

  Eigen::VectorXd residuals;
  int df; // Degrees of freedom

  if (p == 0) {
    // Handle the case when X has zero columns
    if (fit_intercept) {
      // Only intercept model
      double intercept = y.mean();
      residuals = y.array() - intercept;
      df = n - 1;
    } else {
      // Empty model (predicting with zeros)
      residuals = y;
      df = n;
    }
  } else {
    // Normal case with predictors
    auto [ols_intercept, ols_coefs] = fitOls(x, y, fit_intercept);
    residuals = y - x * ols_coefs;

    if (fit_intercept) {
      residuals.array() -= ols_intercept;
    }

    df = n - p - static_cast<int>(fit_intercept);
    assert(df > 0);
  }

  return std::sqrt(residuals.squaredNorm() / df);
}

/**
 * @brief Estimates the regularization parameter alpha for SLOPE regression
 *
 * This function implements an algorithm to estimate an appropriate
 * regularization parameter (alpha) for SLOPE, which is a generalization of the
 * lasso. When n >= p + 30, it directly estimates alpha from OLS residuals.
 * Otherwise, it uses an iterative procedure that alternates between estimating
 * alpha and fitting the SLOPE model.
 *
 * The iterative procedure works by:
 * 1. Starting with an empty set of selected variables
 * 2. Estimating alpha based on the selected variables
 * 3. Fitting a SLOPE model with that alpha
 * 4. Updating the selected variables based on non-zero coefficients
 * 5. Repeating until convergence or maximum iterations reached
 *
 * @tparam MatrixType The type of matrix used to store the design matrix
 * @param x Design matrix with n observations and p predictors
 * @param y Response matrix
 * @param model The SLOPE model object containing method parameters
 * @return A SlopePath object containing the fitted model with estimated alpha
 * @throws std::runtime_error If maximum iterations reached or if too many
 * variables selected
 */
template<typename MatrixType>
SlopePath
estimateAlpha(MatrixType& x, Eigen::MatrixXd& y, Slope& model)
{
  int n = x.rows();
  int p = x.cols();

  int alpha_est_maxit = model.getAlphaEstimationMaxIterations();

  // TODO: Possibly we could just run the path instead. We could use the
  // selected predictors at its end as the selected set and run these in an OLS.
  // There would no longer be any problem with cyclic behavior of this
  // algorithm, but on the other hand the model is likely very unstable.
  std::vector<int> selected;

  Eigen::ArrayXd alpha(1);

  model.setAlphaType("path"); // Otherwise we would be in an infinite loop

  SlopePath result;

  // Estimate the noise level, if possible
  if (n >= p + 30) {
    alpha(0) = estimateNoise(x, y, model.getFitIntercept()) / n;
    result = model.path(x, y, alpha);
  } else {
    for (int it = 0; it < alpha_est_maxit; ++it) {

      MatrixType x_selected = subsetCols(x, selected);

      std::vector<int> selected_prev = selected;
      selected.clear();

      alpha(0) = estimateNoise(x_selected, y, model.getFitIntercept()) / n;

      // TODO: If changes in alpha are small between two steps, then it should
      // be easy to screen, but we are not using that possibility here.
      result = model.path(x, y, alpha);
      auto coefs = result.getCoefs().back();

      for (typename Eigen::SparseMatrix<double>::InnerIterator it(coefs, 0); it;
           ++it) {
        selected.emplace_back(it.row());
      }

      if (selected == selected_prev) {
        return result;
      }

      if (static_cast<int>(selected.size()) >= n + model.getFitIntercept()) {
        throw std::runtime_error(
          "selected >= n - 1 variables, cannot estimate variance");
      }
    }

    slope::WarningLogger::addWarning(
      slope::WarningCode::MAXIT_REACHED,
      "Maximum iterations reached in alpha estimation");
  }

  return result;
}
} // namespace slope
