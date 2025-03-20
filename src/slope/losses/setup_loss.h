/**
 * @file
 * @brief Factory function to create the appropriate loss function based on
 */

#pragma once

#include "loss.h"
#include <memory>

namespace slope {

/**
 * @brief Factory function to create the appropriate loss function based on
 * the distribution family.
 *
 * @details This function creates and returns an loss function object based
 * on the specified statistical distribution family. The supported families are:
 * - "logistic": For binary classification problems (logistic regression)
 * - "poisson": For count data modeling (Poisson regression)
 * - "multinomial": For multi-class classification problems
 * - "quadratic": For continuous response variables (linear regression, default
 *   if unspecified)
 *
 * @param loss A string specifying the loss type ("logistic",
 * "poisson", "multinomial", or "quadratic")
 * @return std::unique_ptr<Loss> A unique pointer to the appropriate
 * loss function object
 */
std::unique_ptr<Loss>
setupLoss(const std::string& loss);

}
