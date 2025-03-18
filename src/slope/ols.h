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
#include <utility>

/**
 * @brief Fits an OLS regression model using dense matrices
 *
 * @param X The design matrix with predictors (features) - dense representation
 * @param y The target/response vector
 * @param fit_intercept Whether to include an intercept term (default: true)
 *
 * @return std::pair containing:
 *         - double: Intercept term (or 0 if fit_intercept is false)
 *         - Eigen::VectorXd: Coefficient vector
 */
std::pair<double, Eigen::VectorXd>
fitOls(const Eigen::MatrixXd& X,
       const Eigen::VectorXd& y,
       bool fit_intercept = true);

/**
 * @brief Fits an OLS regression model using sparse matrices
 *
 * @param X The design matrix with predictors (features) - sparse representation
 * @param y The target/response vector
 * @param fit_intercept Whether to include an intercept term (default: true)
 *
 * @return std::pair containing:
 *         - double: Intercept term (or 0 if fit_intercept is false)
 *         - Eigen::VectorXd: Coefficient vector
 */
std::pair<double, Eigen::VectorXd>
fitOls(const Eigen::SparseMatrix<double>& X,
       const Eigen::VectorXd& y,
       bool fit_intercept = true);
