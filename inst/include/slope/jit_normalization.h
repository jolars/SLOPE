/**
 * @file
 * @brief Enums to control predictor standardization behavior
 */

#pragma once

namespace slope {
/**
 * @brief Enums to control predictor standardization behavior
 */
enum class JitNormalization
{
  None = 0,   ///< No JIT normalization
  Center = 1, ///< Center JIT
  Scale = 2,  ///< Scale JIT
  Both = 3    ///< Both
};

} // namespace slope
