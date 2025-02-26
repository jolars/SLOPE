/**
 * @internal
 * @file
 * @brief Definitions of constants used in libslope
 */

#pragma once

namespace slope {
namespace constants {

constexpr double EPSILON = 1e-10;
constexpr double P_MIN = 1e-9;
constexpr double P_MAX = 1.0 - P_MIN;

} // namespace constants
} // namespace slope
