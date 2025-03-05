/**
 * @file screening.h
 * @brief Screening rules for SLOPE regression optimization
 *
 * Implements feature screening methods to identify active and strong sets
 * of variables, reducing computational complexity in coordinate descent.
 */

#pragma once

#include <Eigen/Core>
#include <vector>

namespace slope {

/**
 * @brief Identifies previously active variables
 *
 * @param beta Current coefficient matrix
 * @return std::vector<int> Indices of variables with non-zero coefficients
 */
std::vector<int>
activeSet(const Eigen::VectorXd& beta);

/**
 * @brief Determines the strong set using sequential strong rules
 *
 * @param gradient_prev Gradient from previous solution
 * @param lambda Current lambda sequence
 * @param lambda_prev Previous lambda sequence
 * @return std::vector<int> Indices of variables in the strong set
 */
std::vector<int>
strongSet(const Eigen::VectorXd& gradient_prev,
          const Eigen::ArrayXd& lambda,
          const Eigen::ArrayXd& lambda_prev);

}
