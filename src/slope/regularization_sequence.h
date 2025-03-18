/**
 * @file
 * @brief Functions for generating regularization sequences for SLOPE.
 * */

#pragma once

#include <Eigen/SparseCore>
#include <string>

namespace slope {

/**
 * Generates a sequence of regularization weights for the sorted L1 norm.
 *
 * @param p The number of lambda values to generate (number of features)
 * @param q The false discovery rate (FDR) level or quantile value (in (0, 1))
 * @param type The type of sequence to generate:
 *            - "bh": Benjamini-Hochberg sequence
 *            - "gaussian": Gaussian sequence
 *            - "oscar": Octagonal Shrinkage and Clustering Algorithm for
 * Regression
 * @param n Number of observations (only used for gaussian type)
 * @param theta1 First parameter for OSCAR weights (default: 1.0)
 * @param theta2 Second parameter for OSCAR weights (default: 1.0)
 * @return Eigen::ArrayXd containing the generated lambda sequence in decreasing
 * order
 */
Eigen::ArrayXd
lambdaSequence(const int p,
               const double q,
               const std::string& type,
               const int n = -1,
               const double theta1 = 1.0,
               const double theta2 = 1.0);

/**
 * Computes a sequence of regularization weights for the SLOPE path.
 *
 * @param alpha_in Alpha sequence, of length zero if automatic.
 * @param path_length Length of path.
 * @param alpha_min_ratio Ratio of minimum to maximum alpha
 * @param alpha_max Value of alpha as which the model is completely sparse
 * @return Eigen::ArrayXd containing the sequence of regularization parameters
 * from strongest (alpha_max) to weakest (alpha_max * alpha_min_ratio)
 */
Eigen::ArrayXd
regularizationPath(const Eigen::ArrayXd& alpha_in,
                   const int path_length,
                   double alpha_min_ratio,
                   const double alpha_max);

} // namespace slope
