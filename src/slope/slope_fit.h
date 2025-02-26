#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>

namespace slope {

/**
 * @brief A class representing the results of SLOPE (Sorted L1 Penalized
 * Estimation) fitting.
 *
 * This class stores the results of a SLOPE regression, including coefficients,
 * intercepts, regularization parameters, and optimization metrics.
 */
class SlopeFit
{
private:
  Eigen::VectorXd intercepts;        ///< Vector of intercept terms
  Eigen::SparseMatrix<double> coefs; ///< Sparse matrix of fitted coefficients
  double alpha;                      ///< Scaling of lambda sequence
  Eigen::ArrayXd lambda; ///< Regularization weights for the sorted L1 norm
  double deviance;       ///< Final model deviance
  double null_deviance;  ///< Null (or intercept-only) model deviance
  std::vector<double>
    primals; ///< History of primal objective values during optimization
  std::vector<double>
    duals; ///< History of dual objective values during optimization
  std::vector<double> time; ///< Time points during optimization
  int passes; ///< Number of passes through the data during optimization

public:
  /// Default constructor
  SlopeFit() = default;

  /**
   * @brief Construct a new Slope Fit object
   *
   * @param intercepts Vector of intercept terms
   * @param coefs Matrix of fitted coefficients
   * @param alpha Mixing parameter between L1 and SLOPE norms
   * @param lambda Sequence of decreasing weights
   * @param deviance Final model deviance
   * @param null_deviance Null model deviance
   * @param primals History of primal objectives
   * @param duals History of dual objectives
   * @param time Vector of optimization timestamps
   * @param passes Number of optimization passes
   */
  SlopeFit(const Eigen::VectorXd& intercepts,
           const Eigen::SparseMatrix<double>& coefs,
           const double alpha,
           const Eigen::ArrayXd& lambda,
           double deviance,
           double null_deviance,
           const std::vector<double>& primals,
           const std::vector<double>& duals,
           const std::vector<double>& time,
           const int passes)
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
   * @brief Gets the intercept terms for this SLOPE fit
   */
  const Eigen::VectorXd& getIntercepts() const { return intercepts; }

  /**
   * @brief Gets the sparse coefficient matrix for this fit
   */
  const Eigen::SparseMatrix<double>& getCoefs() const { return coefs; }

  /**
   * @brief Gets the lambda (regularization) parameter used
   */
  const Eigen::ArrayXd& getLambda() const { return lambda; }

  /**
   * @brief Gets the alpha (mixing) parameter used
   */
  double getAlpha() const { return alpha; }

  /**
   * @brief Gets the model deviance
   */
  double getDeviance() const { return deviance; }

  /**
   * @brief Gets the null model deviance
   */
  double getNullDeviance() const { return null_deviance; }

  /**
   * @brief Gets the sequence of primal objective values during optimization
   */
  const std::vector<double>& getPrimals() const { return primals; }

  /**
   * @brief Gets the sequence of dual objective values during optimization
   */
  const std::vector<double>& getDuals() const { return duals; }

  /**
   * @brief Gets the sequence of computation times during optimization
   */
  const std::vector<double>& getTime() const { return time; }

  /**
   * @brief Gets the total number of optimization iterations
   */
  int getPasses() const { return passes; }

  /**
   * @brief Calculate the deviance ratio (1 - deviance/null_deviance)
   *
   * @return double The deviance ratio, a measure of model fit (higher is
   * better)
   */
  double getDevianceRatios() const { return 1.0 - deviance / null_deviance; }

  /**
   * @brief Calculate the duality gaps during optimization
   *
   * @return std::vector<double> Vector of duality gaps (primal - dual
   * objectives)
   */
  std::vector<double> getGaps() const
  {
    std::vector<double> gaps(primals.size());
    for (size_t i = 0; i < primals.size(); i++) {
      gaps[i] = primals[i] - duals[i];
    }
    return gaps;
  }
};

} // namespace slope
