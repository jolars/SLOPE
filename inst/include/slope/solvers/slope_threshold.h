/**
 * @file
 * @brief The declaration of the slopeThreshold function.
 */

#pragma once

#include "../clusters.h"
#include <Eigen/Core>
#include <tuple>

namespace slope {

/**
 * Calculates slope thresholding for a given input
 *
 * This function calculates the slope threshold for a given x and j, using the
 * provided lambdas and clusters. This is used in the coordinate descent update
 * step.
 *
 * @param x The value of x.
 * @param j The value of j.
 * @param lambda_cumsum Cumulative sum of the lambda sequence.
 * @param clusters The clusters object.
 * @return A tuple containing the slope threshold and the index.
 */
std::tuple<double, int>
slopeThreshold(const double x,
               const int j,
               const Eigen::ArrayXd& lambda_cumsum,
               const Clusters& clusters);

/**
 * Calculates SLOPE thresholding without scaling the lambda sequence.
 *
 * This overload evaluates the thresholding operator in terms of its unscaled
 * linear and quadratic coefficients. It avoids materializing a scaled copy of
 * the cumulative lambda sequence for every coordinate update.
 *
 * @param gamma The linear coefficient of the one-dimensional problem.
 * @param omega The positive quadratic coefficient.
 * @param j The index of the cluster being updated.
 * @param lambda_cumsum Cumulative sum of the unscaled lambda sequence.
 * @param clusters The clusters object.
 * @return A tuple containing the slope threshold and the index.
 */
std::tuple<double, int>
slopeThreshold(const double gamma,
               const double omega,
               const int j,
               const Eigen::ArrayXd& lambda_cumsum,
               const Clusters& clusters);

}
