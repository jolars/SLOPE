/**
 * @file
 * @brief Ordinary Least Squares (OLS) regression functionality
 *
 * This header provides functions for fitting Ordinary Least Squares regression
 * models using both dense and sparse matrix representations.
 */

#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <utility>

namespace slope {
namespace detail {

/**
 * Fits ordinary least squares (OLS) regression using dense matrices.
 *
 * This function performs OLS regression using QR decomposition with column
 * pivoting for numerical stability. Optionally includes an intercept term.
 *
 * @param X Feature matrix (dense)
 * @param y Response vector
 * @param fit_intercept Whether to include an intercept term (default: true)
 *
 * @return std::pair<double, Eigen::VectorXd> containing:
 *         - first: Intercept value (0.0 if fit_intercept is false)
 *         - second: Coefficient vector for predictors
 */
template<typename T>
std::pair<double, Eigen::VectorXd>
fitOls(const Eigen::MatrixBase<T>& X,
       const Eigen::VectorXd& y,
       bool fit_intercept = true)
{
  Eigen::MatrixXd x_mod = X;

  if (fit_intercept) {
    // Add column of ones for intercept
    Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(x_mod.rows(), 1);
    x_mod.conservativeResize(Eigen::NoChange, x_mod.cols() + 1);
    x_mod.rightCols(1) = ones;
  }

  // Solve with column pivoting
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(x_mod);
  Eigen::VectorXd all_coefs = qr.solve(y);

  double intercept = 0.0;
  Eigen::VectorXd coeffs;

  if (fit_intercept) {
    intercept = all_coefs(all_coefs.size() - 1);
    coeffs = all_coefs.head(all_coefs.size() - 1);
  } else {
    coeffs = all_coefs;
  }

  return { intercept, coeffs };
}

/**
 * Fits ordinary least squares (OLS) regression using sparse matrices.
 *
 * This function performs OLS regression for sparse feature matrices using
 * sparse matrix operations. Optionally includes an intercept term by
 * augmenting the sparse matrix with a dense column of ones.
 *
 * @param X Feature matrix (sparse)
 * @param y Response vector
 * @param fit_intercept Whether to include an intercept term (default: true)
 *
 * @return std::pair<double, Eigen::VectorXd> containing:
 *         - first: Intercept value (0.0 if fit_intercept is false)
 *         - second: Coefficient vector for predictors
 */
template<typename T>
std::pair<double, Eigen::VectorXd>
fitOls(const Eigen::SparseMatrixBase<T>& X,
       const Eigen::VectorXd& y,
       bool fit_intercept = true)
{
  // TODO: Investigate if we can avoid this copy.
  Eigen::SparseMatrix<double> x_mod = X;

  if (fit_intercept) {
    // Construct column of ones for intercept
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x_mod.rows());
    Eigen::SparseMatrix<double> intercept_col(x_mod.rows(), 1);
    intercept_col.reserve(x_mod.rows());

    for (int i = 0; i < x_mod.rows(); ++i) {
      intercept_col.insert(i, 0) = ones(i);
    }

    // Concatenate X and intercept column
    Eigen::SparseMatrix<double> temp(x_mod.rows(), x_mod.cols() + 1);
    temp.leftCols(x_mod.cols()) = x_mod;
    temp.rightCols(1) = intercept_col;
    x_mod = temp;
  }

  // Solve with sparse QR
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
    solver;
  solver.compute(x_mod);
  Eigen::VectorXd all_cofs = solver.solve(y);

  double intercept = 0.0;
  Eigen::VectorXd coeffs;

  if (fit_intercept) {
    intercept = all_cofs(all_cofs.size() - 1);
    coeffs = all_cofs.head(all_cofs.size() - 1);
  } else {
    coeffs = all_cofs;
  }

  return { intercept, coeffs };
}

} // namespace detail
} // namespace slope
