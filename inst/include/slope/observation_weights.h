#pragma once

#include <Eigen/Core>
#include <stdexcept>

namespace slope::detail {

inline bool
uniformObservationWeights(const Eigen::VectorXd& weights)
{
  return weights.size() == 0 ||
         (weights(0) > 0.0 && (weights.array() == weights(0)).all());
}

inline void
validateObservationWeights(const Eigen::VectorXd& weights, Eigen::Index n)
{
  if (n <= 0) {
    throw std::invalid_argument("Observation weights require at least one row");
  }
  if (weights.size() == 0) {
    return;
  }
  if (weights.size() != n) {
    throw std::invalid_argument(
      "Observation weights must have one entry per row");
  }
  if (!weights.allFinite() || (weights.array() < 0.0).any()) {
    throw std::invalid_argument(
      "Observation weights must be finite and nonnegative");
  }
  if (weights.maxCoeff() == 0.0) {
    throw std::invalid_argument(
      "Observation weights must have positive total weight");
  }
}

/** Validate and normalize row weights to have mean one. */
inline Eigen::VectorXd
normalizeObservationWeights(const Eigen::VectorXd& weights, Eigen::Index n)
{
  validateObservationWeights(weights, n);
  if (uniformObservationWeights(weights)) {
    return Eigen::VectorXd::Ones(n);
  }
  // Scaling first avoids overflow when finite weights have a very large sum.
  Eigen::VectorXd normalized = weights / weights.maxCoeff();
  normalized *= static_cast<double>(n) / normalized.sum();
  return normalized;
}

inline void
validateBinomialClassWeights(const Eigen::MatrixXd& y,
                             const Eigen::VectorXd& weights)
{
  for (Eigen::Index j = 0; j < y.cols(); ++j) {
    const auto positive_weight = weights.array() > 0.0;
    if (!(positive_weight && (y.col(j).array() == 0.0)).any() ||
        !(positive_weight && (y.col(j).array() == 1.0)).any()) {
      throw std::invalid_argument(
        "Binomial fits with an intercept require positive total weight in "
        "both classes");
    }
  }
}

inline void
validateMultinomialClassWeights(const Eigen::MatrixXd& y,
                                const Eigen::VectorXd& weights)
{
  const auto positive_weight = weights.array() > 0.0;
  const bool reference_present =
    (positive_weight && (y.rowwise().sum().array() == 0.0)).any();
  if (!reference_present || ((weights.transpose() * y).array() <= 0.0).any()) {
    throw std::invalid_argument(
      "Multinomial fits with an intercept require positive total weight in "
      "every class");
  }
}

inline Eigen::Index
positiveWeightCount(const Eigen::VectorXd& weights)
{
  return (weights.array() > 0.0).count();
}

inline void
maskZeroWeightRows(Eigen::MatrixXd& values, const Eigen::VectorXd& weights)
{
  for (Eigen::Index i = 0; i < weights.size(); ++i) {
    if (weights(i) == 0.0) {
      values.row(i).setZero();
    }
  }
}

inline Eigen::RowVectorXd
weightedMean(const Eigen::MatrixXd& values, const Eigen::VectorXd& weights)
{
  Eigen::RowVectorXd mean = Eigen::RowVectorXd::Zero(values.cols());
  const double total = weights.sum();
  for (Eigen::Index i = 0; i < values.rows(); ++i) {
    if (weights(i) > 0.0) {
      mean += (weights(i) / total) * values.row(i);
    }
  }
  return mean;
}

} // namespace slope::detail
