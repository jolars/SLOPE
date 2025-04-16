/**
 * @file
 * @brief SLOPE (Sorted L-One Penalized Estimation) fitting results
 */

#pragma once

#include "clusters.h"
#include "losses/setup_loss.h"
#include "normalize.h"
#include "utils.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <memory>

namespace slope {

/**
 * @brief A class representing the results of SLOPE (Sorted L1 Penalized
 * Estimation) fitting.
 *
 * This class stores the results of a SLOPE regression, including coefficients,
 * intercepts, clusters (if required), regularization parameters, and
 * optimization metrics.
 */
class SlopeFit
{
private:
  Eigen::VectorXd intercepts;        ///< Vector of intercept terms
  Eigen::SparseMatrix<double> coefs; ///< Sparse matrix of fitted coefficients
  Clusters clusters;                 ///< Clusters
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
  std::string centering_type; ///< Type of centering
  std::string scaling_type;   ///< Type of scaling
  std::string loss_type;      ///< Loss type

public:
  SlopeFit() = default;

  /**
   * @brief Construct a new Slope Fit object
   *
   * @param intercepts Vector of intercept terms
   * @param coefs Matrix of fitted coefficients
   * @param clusters Clusters of coefficients
   * @param alpha Mixing parameter between L1 and SLOPE norms
   * @param lambda Sequence of decreasing weights
   * @param deviance Final model deviance
   * @param null_deviance Null model deviance
   * @param primals History of primal objectives
   * @param duals History of dual objectives
   * @param time Vector of optimization timestamps
   * @param passes Number of optimization passes
   * @param centering_type Type of centering for the design matrix
   * @param scaling_type Type of scaling for the design matrix
   *
   */
  SlopeFit(const Eigen::VectorXd& intercepts,
           const Eigen::SparseMatrix<double>& coefs,
           const Clusters& clusters,
           const double alpha,
           const Eigen::ArrayXd& lambda,
           const double deviance,
           const double null_deviance,
           const std::vector<double>& primals,
           const std::vector<double>& duals,
           const std::vector<double>& time,
           const int passes,
           const std::string& centering_type,
           const std::string& scaling_type)
    : intercepts{ intercepts }
    , coefs{ coefs }
    , clusters{ clusters }
    , alpha{ alpha }
    , lambda{ lambda }
    , deviance{ deviance }
    , null_deviance{ null_deviance }
    , primals{ primals }
    , duals{ duals }
    , time{ time }
    , passes{ passes }
    , centering_type{ centering_type }
    , scaling_type{ scaling_type }
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
   * @brief Gets the clusters
   */
  const Clusters& getClusters() const { return clusters; }

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
   * @brief Gets the model's loss type
   */
  const std::string& getLossType() const { return loss_type; }

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
  double getDevianceRatio() const { return 1.0 - deviance / null_deviance; }

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

  /**
   * @brief Predict the response for a given input matrix
   * @tparam T Type of input matrix (dense or sparse)
   * @param x Input matrix of features
   * @param type Type of prediction to return ("response" or "linear")
   * @return Matrix of predicted responses
   * @see Loss
   */
  template<typename T>
  Eigen::MatrixXd predict(T& x, const std::string& type = "response") const
  {
    validateOption(type, { "response", "linear" }, "type");

    int n = x.rows();
    int m = coefs.cols();

    std::unique_ptr<Loss> loss = setupLoss(this->loss_type);

    Eigen::VectorXd x_centers;
    Eigen::VectorXd x_scales;
    normalize(
      x, x_centers, x_scales, this->centering_type, this->scaling_type, false);

    Eigen::MatrixXd eta = Eigen::MatrixXd::Zero(n, m);

    for (int k = 0; k < m; ++k) {
      for (typename Eigen::SparseMatrix<double>::InnerIterator it(coefs, k); it;
           ++it) {
        int j = it.row();
        eta.col(k) += x.col(j) * (it.value() / x_scales(j));
        eta.col(k).array() -= it.value() * x_centers(j) / x_scales(j);
      }

      eta.col(k).array() += intercepts(k);
    }

    if (type == "linear") {
      return eta;
    }

    // Return predictions
    return loss->predict(eta);
  }
};

} // namespace slope
