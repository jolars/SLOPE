/**
 * @file
 * @brief An implementation of the quantile function for the standard normal
 * distribution
 */

#include <cmath>

namespace slope {

/**
 * @brief Computes the quantile of a standard normal distribution using the
 * Beasley-Springer-Moro algorithm.
 *
 * @param p Probability value
 * @return Quantile
 */
double
normalQuantile(const double p);

} // namespace slope
