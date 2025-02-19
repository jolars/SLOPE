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

#include <Eigen/Dense>
#include <Eigen/SparseCore>

namespace slope {

/**
 * @class SlopePath
 * @brief Container class for SLOPE regression solution paths
 *
 * Stores and provides access to:
 * - Model coefficients (intercepts and sparse coefficient matrices)
 * - Regularization parameters (alpha and lambda sequences)
 * - Model fit statistics (deviance, null deviance)
 * - Optimization metrics (primal/dual objectives, computation time, iteration
 * counts)
 */
class SlopePath
{
private:
  std::vector<Eigen::VectorXd> intercepts;
  std::vector<Eigen::SparseMatrix<double>> coefs;
  Eigen::ArrayXd alpha;
  Eigen::ArrayXd lambda;
  std::vector<double> deviance;
  double null_deviance;
  std::vector<std::vector<double>> primals;
  std::vector<std::vector<double>> duals;
  std::vector<std::vector<double>> time;
  std::vector<int> passes;

public:
  /**
   * @brief Constructs an empty SlopePath object
   */
  SlopePath() = default;

  /**
   * @brief Constructs a SlopePath object containing SLOPE regression solution
   * path data
   *
   * @param intercepts Vector of intercept terms for each solution
   * @param coefs Vector of sparse coefficient matrices for each solution
   * @param alpha Vector of mixing parameters (elastic net mixing parameter)
   * @param lambda Vector of regularization parameters
   * @param deviance Vector of model deviances for each solution
   * @param null_deviance Null model deviance (deviance of intercept-only model)
   * @param primals Vector of primal objective values during optimization
   * @param duals Vector of dual objective values during optimization
   * @param time Vector of computation times for each solution
   * @param passes Vector of iteration counts for each solution
   *
   * @note All vectors should have consistent lengths corresponding to the
   * number of solutions in the path
   */
  SlopePath(const std::vector<Eigen::VectorXd>& intercepts,
            const std::vector<Eigen::SparseMatrix<double>>& coefs,
            const Eigen::ArrayXd& alpha,
            const Eigen::ArrayXd& lambda,
            const std::vector<double>& deviance,
            double null_deviance,
            const std::vector<std::vector<double>>& primals,
            const std::vector<std::vector<double>>& duals,
            const std::vector<std::vector<double>>& time,
            const std::vector<int>& passes)
    : intercepts{ intercepts }
    , coefs{ coefs }
    , alpha{ alpha }
    , lambda{ lambda }
    , deviance{ deviance }
    , null_deviance{ null_deviance }
    , primals{ primals }
    , duals{ duals }
    , time{ time }
    , passes{ passes }
  {
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
  const std::vector<Eigen::VectorXd>& getIntercepts() const
  {
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
  const std::vector<Eigen::SparseMatrix<double>>& getCoefs() const
  {
    return coefs;
  }

  /**
   * @brief Gets the sparse coefficient matrix for a specific solution in the
   * path
   *
   * @param i Index of the solution
   * @return const Eigen::SparseMatrix<double>& Reference to the coefficient
   * matrix at index i
   * @throws std::runtime_error if index is out of bounds
   */
  const Eigen::SparseMatrix<double>& getCoefs(const std::size_t i) const
  {
    assert(i < coefs.size() && "Index out of bounds");
    return coefs[i];
  }

  /**
   * @brief Gets the alpha parameter sequence
   */
  const Eigen::ArrayXd& getAlpha() const { return alpha; }

  /**
   * @brief Gets the lambda (regularization) weights
   */
  const Eigen::ArrayXd& getLambda() const { return lambda; }

  /**
   * @brief Gets the deviance values for each solution
   */
  const std::vector<double>& getDeviance() const { return deviance; }

  /**
   * @brief Gets the null model deviance
   */
  double getNullDeviance() const { return null_deviance; }

  /**
   * @brief Gets the primal objective values during optimization
   */
  const std::vector<std::vector<double>>& getPrimals() const { return primals; }

  /**
   * @brief Gets the dual objective values during optimization
   */
  const std::vector<std::vector<double>>& getDuals() const { return duals; }

  /**
   * @brief Gets the computation times for each solution
   */
  const std::vector<std::vector<double>>& getTime() const { return time; }

  /**
   * @brief Gets the number of iterations for each solution
   */
  const std::vector<int>& getPasses() const { return passes; }

  /**
   * @brief Computes the deviance ratio (explained deviance) for each solution
   * @return std::vector<double> Vector of deviance ratios (1 -
   * deviance/null_deviance)
   */
  const std::vector<double> getDevianceRatios() const
  {
    std::vector<double> ratios(deviance.size());

    for (size_t i = 0; i < deviance.size(); i++) {
      ratios[i] = 1.0 - deviance[i] / null_deviance;
    }

    return ratios;
  }

  /**
   * @brief Computes the duality gaps (primal - dual objectives) for each
   * solution
   * @return std::vector<std::vector<double>> Vector of duality gap sequences
   */
  const std::vector<std::vector<double>> getGaps() const
  {
    std::vector<std::vector<double>> gaps(primals.size());

    for (size_t i = 0; i < primals.size(); i++) {
      gaps[i].resize(primals[i].size());

      for (size_t j = 0; j < primals[i].size(); j++) {
        gaps[i][j] = primals[i][j] - duals[i][j];
      }
    }

    return gaps;
  }
};

} // namespace slope
