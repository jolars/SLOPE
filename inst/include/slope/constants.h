/**
 * @file
 * @brief Definitions of constants used in libslope
 *
 * This file contains fundamental constants used throughout the slope library.
 * These constants include numerical limits, probability bounds, and special
 * values used for computation and validation.
 */

#pragma once

#include <limits>

namespace slope {

/**
 * @brief Namespace containing constants used in the slope library
 *
 * This namespace contains various constants used in the slope library, such as
 * machine epsilon, minimum and maximum probabilities, and positive and negative
 * infinity.
 */
namespace constants {

/// @brief Small value used for floating-point comparisons to handle precision
/// issues
constexpr double EPSILON = 1e-10;

/// @brief Minimum allowed probability value to avoid numerical underflow
constexpr double P_MIN = 1e-9;

/// @brief Maximum allowed probability value to avoid numerical issues near 1.0
constexpr double P_MAX = 1.0 - P_MIN;

/// @brief Representation of positive infinity using maximum double value
constexpr double POS_INF = std::numeric_limits<double>::max();

/// @brief Representation of negative infinity using lowest double value
constexpr double NEG_INF = std::numeric_limits<double>::lowest();

/// @brief Maximum allowed exponent
constexpr double MAX_EXP = 250;

/// @brief Maximum allowed divisor
constexpr double MAX_DIV = 1e-12;

} // namespace constants
} // namespace slope
