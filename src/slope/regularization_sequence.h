/**
 * @file
 * @brief Functions for generating regularization sequences for SLOPE.
 * */

#pragma once

#include "sorted_l1_norm.h"
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
 * @param gradient The gradient
 * @param penalty Penalty object.
 * @param lambda Regularization weights.
 * @param n Number of observations.
 * @param path_length Length of path.
 * @param alpha_min_ratio Ratio of minimum to maximum alpha
 * @return Eigen::ArrayXd containing the sequence of regularization parameters
 * from strongest (alpha_max) to weakest (alpha_max * alpha_min_ratio)
 */
std::tuple<Eigen::ArrayXd, double, int>
regularizationPath(const Eigen::ArrayXd& alpha_in,
                   const Eigen::VectorXd& gradient,
                   const SortedL1Norm& penalty,
                   const Eigen::ArrayXd& lambda,
                   const int n,
                   const int path_length,
                   double alpha_min_ratio);

} // namespace slope
