/**
 * @file
 * @brief Eigen compatibility layer for version differences
 */

#pragma once

#include <Eigen/Core>

namespace slope {

#if EIGEN_VERSION_AT_LEAST(5, 0, 0)
using Eigen::placeholders::all;
using Eigen::placeholders::last;
#else
using Eigen::all;
using Eigen::last;
#endif

} // namespace slope
