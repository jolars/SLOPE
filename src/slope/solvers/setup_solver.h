/**
 * @file
 * @brief Factory function to create and configure a SLOPE solver
 */

#include "solver.h"
#include <memory>
#include <string>

namespace slope {

/**
 * @brief Factory function to create and configure a SLOPE solver
 *
 * @details Creates a solver object based on the specified type and parameters.
 * The solver implements the Sorted L1 Penalized Estimation (SLOPE) algorithm
 * with various configurations possible.
 *
 * @param solver_type Type of solver to use (e.g., "pgd", "admm")
 * @param loss Loss type
 * @param jit_normalization Type of JIT normalization
 * @param intercept Whether to fit an intercept term
 * @param update_clusters Whether to update cluster assignments during
 * optimization (Hybrid solver)
 * @param cd_iterations Frequency of proximal gradient descent updates (Hybrid
 * solver)
 *
 * @return std::unique_ptr<SolverBase> A unique pointer to the
 * configured solver
 * @see Loss
 * @see JitNormalization
 * @see SolverBase
 * @see PGD
 * @see Hybrid
 */
std::unique_ptr<SolverBase>
setupSolver(const std::string& solver_type,
            const std::string& loss,
            JitNormalization jit_normalization,
            bool intercept,
            bool update_clusters,
            int cd_iterations);

} // namespace slope
