/**
 * @file
 * @brief Karush-Kuhn-Tucker (KKT) optimality condition checking for SLOPE
 * regression
 */

#pragma once

#include <Eigen/Core>

namespace slope {

/**
 * @typedef ArrayXb
 * @brief Dynamic-size column vector of boolean values
 * Wrapper around Eigen::Array<bool, Eigen::Dynamic, 1>
 */
using ArrayXb = Eigen::Array<bool, Eigen::Dynamic, 1>;

/**
 * @typedef ArrayXXb
 * @brief Dynamic-size matrix of boolean values
 * Wrapper around Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>
 */
using ArrayXXb = Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * @brief Checks KKT conditions for SLOPE optimization
 *
 * @param gradient The gradient of the loss function
 * @param beta The current coefficients
 * @param lambda Vector of regularization parameters
 * @param strong_set Vector of indices in the strong set
 * @return std::vector<int> Indices where KKT conditions are violated
 *
 * Verifies if the current solution satisfies the KKT optimality conditions
 * for the SLOPE optimization problem. Returns indices where violations occur.
 */
std::vector<int>
kktCheck(const Eigen::VectorXd& gradient,
         const Eigen::VectorXd& beta,
         const Eigen::ArrayXd& lambda,
         const std::vector<int>& strong_set);

} // namespace slope
