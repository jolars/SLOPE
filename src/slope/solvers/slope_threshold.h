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
 * @param lambdas The array of lambdas.
 * @param clusters The clusters object.
 * @return A tuple containing the slope threshold and the index.
 */
std::tuple<double, int>
slopeThreshold(const double x,
               const int j,
               const Eigen::ArrayXd lambdas,
               const Clusters& clusters);

}
