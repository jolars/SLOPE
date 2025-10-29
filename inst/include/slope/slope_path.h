/**
 * @file
 * @brief Defines the SlopePath class for storing and accessing SLOPE regression
 * solution paths
 *
 * This header file provides the SlopePath class which stores the results of
 * SLOPE (Sorted L-One Penalized Estimation) regression fits across multiple
 * regularization parameters (lambda) and mixing parameters (alpha).
 */

#pragma once

#include "slope_fit.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>

namespace slope {

/**
 * @class SlopePath
 * @brief Container class for SLOPE regression solution paths
 *
 * Stores a path of SlopeFit objects and provides convenience access to:
 * - Model coefficients (intercepts and sparse coefficient matrices)
 * - Regularization parameters (alpha and lambda sequences)
 * - Model fit statistics (deviance, null deviance)
 * - Optimization metrics (primal/dual objectives, computation time, iteration
 * counts)
 */
class SlopePath
{
private:
  std::vector<SlopeFit> fits;

public:
  /**
   * @brief Constructs an empty SlopePath object
   */
  SlopePath() = default;

  /**
   * @brief Constructs a SlopePath object containing SLOPE regression solution
   * path data
   *
   * @param fits Vector of SlopeFit objects for each solution in the path
   */
  SlopePath(const std::vector<SlopeFit>& fits)
    : fits{ fits }
  {
  }

  /**
   * @brief Get one of the coefficients along the path
   *
   * @param step the
   * @return The fit at step `step` of the path. A SlopeFit object.
   */
  const SlopeFit& operator()(const size_t step) const
  {
    assert(step < fits.size());

    return fits[step];
  }

  /**
   * @brief Returns the vector of intercept terms for each solution in the path
   *
   * @return const std::vector<Eigen::VectorXd>& Reference to the vector of
   * intercepts
   *
   * Each element in the returned vector corresponds to the intercept term(s)
   * for a particular solution in the regularization path.
   */
  std::vector<Eigen::VectorXd> getIntercepts(
    const bool original_scale = true) const
  {
    std::vector<Eigen::VectorXd> intercepts;

    for (const auto& fit : fits) {
      intercepts.emplace_back(fit.getIntercepts(original_scale));
    }

    return intercepts;
  }

  /**
   * @brief Returns the vector of coefficient matrices for each solution in the
   * path
   *
   * @return const std::vector<Eigen::SparseMatrix<double>>& Reference to the
   * vector of sparse coefficient matrices
   *
   * Each element in the returned vector is a sparse matrix containing the model
   * coefficients for a particular solution in the regularization path.
   */
  std::vector<Eigen::SparseMatrix<double>> getCoefs(
    const bool original_scale = true) const
  {
    std::vector<Eigen::SparseMatrix<double>> coefs;

    for (const auto& fit : fits) {
      coefs.emplace_back(fit.getCoefs(original_scale));
    }

    return coefs;
  }

  /**
   * @brief Returns the clusters for each solution in the path.
   *
   * @return Vector of clusters
   */
  std::vector<Clusters> getClusters() const
  {
    std::vector<Clusters> clusters;

    for (const auto& fit : fits) {
      clusters.emplace_back(fit.getClusters());
    }

    return clusters;
  }

  /**
   * @brief Gets the alpha parameter sequence
   */
  Eigen::ArrayXd getAlpha() const
  {

    Eigen::ArrayXd alphas(fits.size());

    for (size_t i = 0; i < fits.size(); i++) {
      alphas(i) = fits[i].getAlpha();
    }

    return alphas;
  }

  /**
   * @brief Gets the lambda (regularization) weights
   */
  const Eigen::ArrayXd& getLambda() const { return fits.front().getLambda(); }

  /**
   * @brief Gets the deviance values for each solution
   */
  std::vector<double> getDeviance() const
  {
    std::vector<double> deviances;

    for (const auto& fit : fits) {
      deviances.emplace_back(fit.getDeviance());
    }

    return deviances;
  }

  /**
   * @brief Gets the null model deviance
   */
  double getNullDeviance() const { return fits.front().getNullDeviance(); }

  /**
   * @brief Gets the primal objective values during optimization
   */
  std::vector<std::vector<double>> getPrimals() const
  {
    std::vector<std::vector<double>> primals;

    for (const auto& fit : fits) {
      primals.emplace_back(fit.getPrimals());
    }

    return primals;
  }

  /**
   * @brief Gets the dual objective values during optimization
   */
  std::vector<std::vector<double>> getDuals() const
  {
    std::vector<std::vector<double>> duals;

    for (const auto& fit : fits) {
      duals.emplace_back(fit.getDuals());
    }

    return duals;
  }

  /**
   * @brief Gets the computation times for each solution
   */
  std::vector<std::vector<double>> getTime() const
  {
    std::vector<std::vector<double>> time;

    for (const auto& fit : fits) {
      time.emplace_back(fit.getTime());
    }

    return time;
  }

  /**
   * @brief Gets the number of iterations for each solution
   */
  std::vector<int> getPasses() const
  {
    std::vector<int> passes;

    for (const auto& fit : fits) {
      passes.emplace_back(fit.getPasses());
    }

    return passes;
  }

  /**
   * @brief Computes the deviance ratio (explained deviance) for each solution
   * @return std::vector<double> Vector of deviance ratios (1 -
   * deviance/null_deviance)
   */
  const std::vector<double> getDevianceRatios() const
  {
    std::vector<double> ratios;

    for (const auto& fit : fits) {
      ratios.emplace_back(fit.getDevianceRatio());
    }

    return ratios;
  }

  /**
   * @brief Computes the duality gaps (primal - dual objectives) for each
   * solution
   * @return std::vector<std::vector<double>> Vector of duality gap sequences
   */
  std::vector<std::vector<double>> getGaps() const
  {
    std::vector<std::vector<double>> gaps;

    for (const auto& fit : fits) {
      gaps.emplace_back(fit.getGaps());
    }

    return gaps;
  }

  /**
   * @brief Gets the number of solutions in the path
   * @return Size of the path (number of SlopeFit objects)
   */
  std::size_t size() const { return fits.size(); }
};

} // namespace slope
