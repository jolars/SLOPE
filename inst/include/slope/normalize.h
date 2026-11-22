/** @file
 * @brief Functions to normalize the design matrix and rescale coefficients
 * in case the design was normalized
 */

#pragma once

#include "jit_normalization.h"
#include "math.h"
#include "observation_weights.h"
#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <limits>

namespace slope {

namespace detail {

class WeightRangeSums
{
public:
  explicit WeightRangeSums(const Eigen::VectorXd& weights)
    : weights(weights)
  {
  }

  double sum(Eigen::Index begin, Eigen::Index end);

private:
  const Eigen::VectorXd& weights;
  Eigen::VectorXd sums;
};

template<typename T, typename F>
void
forEachWeightedColumnValue(const Eigen::MatrixBase<T>& x,
                           const Eigen::VectorXd& weights,
                           const int column,
                           const Eigen::Index,
                           WeightRangeSums*,
                           F visit)
{
  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    if (weights(i) > 0.0) {
      visit(x(i, column), weights(i));
    }
  }
}

template<typename T, typename F>
void
forEachWeightedColumnValue(const Eigen::SparseMatrixBase<T>& x,
                           const Eigen::VectorXd& weights,
                           const int column,
                           const Eigen::Index support,
                           WeightRangeSums* range_sums,
                           F visit)
{
  Eigen::Index stored_rows = 0;
  for (typename T::InnerIterator it(x.derived(), column); it; ++it) {
    if (weights(it.row()) > 0.0) {
      visit(it.value(), weights(it.row()));
      ++stored_rows;
    }
  }

  // Implicit zeros matter only when they occur on positive-weight rows.
  if (stored_rows < support) {
    double implicit_weight = 0.0;
    if (range_sums) {
      Eigen::Index next_row = 0;
      for (typename T::InnerIterator it(x.derived(), column); it; ++it) {
        implicit_weight += range_sums->sum(next_row, it.row());
        next_row = it.row() + 1;
      }
      implicit_weight += range_sums->sum(next_row, weights.size());
    }
    visit(0.0, implicit_weight);
  }
}

template<typename T>
Eigen::VectorXd
weightedFeatureCenters(const T& x,
                       const std::string& type,
                       const Eigen::VectorXd& weights)
{
  Eigen::VectorXd centers(x.cols());
  const double total = weights.sum();
  const Eigen::Index support = positiveWeightCount(weights);
  for (Eigen::Index j = 0; j < x.cols(); ++j) {
    double value =
      type == "min" ? std::numeric_limits<double>::infinity() : 0.0;
    forEachWeightedColumnValue(
      x, weights, j, support, nullptr, [&](double entry, double weight) {
        if (type == "min") {
          value = std::min(value, entry);
        } else {
          value += (weight / total) * entry;
        }
      });
    centers(j) = value;
  }
  return centers;
}

template<typename T>
Eigen::VectorXd
weightedFeatureScales(const T& x,
                      const std::string& type,
                      const Eigen::VectorXd& weights)
{
  Eigen::VectorXd scales(x.cols());
  Eigen::VectorXd means;
  if (type == "sd") {
    means = weightedFeatureCenters(x, "mean", weights);
  }
  const double total = weights.sum();
  const Eigen::Index support = positiveWeightCount(weights);
  WeightRangeSums range_sums(weights);
  for (Eigen::Index j = 0; j < x.cols(); ++j) {
    double value = 0.0;
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    forEachWeightedColumnValue(
      x,
      weights,
      j,
      support,
      type == "sd" ? &range_sums : nullptr,
      [&](double entry, double weight) {
        if (type == "range") {
          minimum = std::min(minimum, entry);
          maximum = std::max(maximum, entry);
        } else if (type == "max_abs") {
          value = std::max(value, std::abs(entry));
        } else if (type == "l1") {
          value += (weight / total) * std::abs(entry);
        } else {
          double centered = entry;
          if (type == "sd") {
            minimum = std::min(minimum, entry);
            maximum = std::max(maximum, entry);
            centered -= means(j);
          }
          value = std::hypot(value, std::sqrt(weight / total) * centered);
        }
      });
    if (type == "range") {
      value = maximum - minimum;
    } else if (type == "sd" && minimum == maximum) {
      // Rounding in the mean must not give constant columns a nonzero scale.
      value = 0.0;
    } else if (type == "l1") {
      value *= support;
    } else if (type == "l2") {
      value *= std::sqrt(support);
    }
    scales(j) = value;
  }
  return scales;
}

} // namespace detail

/**
 * Compute centers.
 *
 * Supported types are "none", "manual", "mean", and "min". With observation
 * weights, means use normalized weights, and minima ignore zero-weight rows.
 *
 * @tparam T The type of the input matrix.
 * @param x_centers A vector where the computed or provided centers will be
 * stored.
 * @param x The input matrix.
 * @param type A string specifying the normalization type ("none", "manual", or
 * "standardization").
 * @param observation_weights Optional nonnegative weights, one per row.
 *
 * @throws std::invalid_argument if the provided manual centers or scales have
 * invalid dimensions or contain non-finite values.
 */
template<typename T>
void
computeCenters(Eigen::VectorXd& x_centers,
               const T& x,
               const std::string& type,
               const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
{
  int p = x.cols();
  Eigen::VectorXd weights;
  if (observation_weights.size() > 0) {
    weights =
      detail::normalizeObservationWeights(observation_weights, x.rows());
    if (detail::uniformObservationWeights(weights)) {
      weights.resize(0);
    }
  }

  if (type == "manual") {
    if (x_centers.size() != p) {
      throw std::invalid_argument("Invalid dimensions in centers");
    }

    if (!x_centers.allFinite()) {
      throw std::invalid_argument("Centers must be finite");
    }

  } else if (weights.size() > 0 && (type == "mean" || type == "min")) {
    x_centers = detail::weightedFeatureCenters(x, type, weights);
  } else if (type == "mean") {
    x_centers = means(x);
  } else if (type == "min") {
    x_centers = mins(x);
  } else if (type != "none") {
    throw std::invalid_argument("Invalid centering type");
  }
}

/**
 * Compute scales
 *
 * Supported types are "none", "manual", "sd", "l1", "l2", "max_abs", and
 * "range". Weighted standard deviations use the population denominator.
 * Weighted L1 and L2 norms normalize weights to the positive-weight row count
 * to preserve the unweighted norm convention. Extrema ignore zero-weight rows.
 *
 * @tparam T The type of the input matrix.
 * @param x The input matrix.
 * @param x_scales A vector where the computed or provided scales will be
 * stored.
 * @param type A string specifying the normalization type ("none", "manual", or
 * "standardization").
 * @param observation_weights Optional nonnegative weights, one per row.
 *
 * @throws std::invalid_argument if the provided manual centers or scales have
 * invalid dimensions or contain non-finite values.
 */
template<typename T>
void
computeScales(Eigen::VectorXd& x_scales,
              const T& x,
              const std::string& type,
              const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
{
  int p = x.cols();
  Eigen::VectorXd weights;
  if (observation_weights.size() > 0) {
    weights =
      detail::normalizeObservationWeights(observation_weights, x.rows());
    if (detail::uniformObservationWeights(weights)) {
      weights.resize(0);
    }
  }

  if (type == "manual") {
    if (x_scales.size() != p) {
      throw std::invalid_argument("Invalid dimensions in scales");
    }
    if (!x_scales.allFinite()) {
      throw std::invalid_argument("Scales must be finite");
    }
  } else if (weights.size() > 0 &&
             (type == "sd" || type == "l1" || type == "l2" ||
              type == "max_abs" || type == "range")) {
    x_scales = detail::weightedFeatureScales(x, type, weights);
  } else if (type == "sd") {
    x_scales = stdDevs(x);
  } else if (type == "l1") {
    x_scales = l1Norms(x);
  } else if (type == "l2") {
    x_scales = l2Norms(x);
  } else if (type == "max_abs") {
    x_scales = maxAbs(x);
  } else if (type == "range") {
    x_scales = ranges(x);
  } else if (type != "none") {
    throw std::invalid_argument("Invalid scaling type");
  }
}

/**
 * Normalize a dense matrix by centering and scaling.
 *
 * The function computes column centers and scaling factors based on the
 * specified normalization type ("none", "manual", or "standardization"). If
 * modify_x is true, the normalization is applied directly to the input matrix.
 *
 * @param x The dense input matrix to be normalized.
 * @param x_centers A vector that will hold the column centers. It will be
 * resized to match the number of columns.
 * @param x_scales  A vector that will hold the column scaling factors. It will
 * be resized to match the number of columns.
 * @param centering_type A string specifying the normalization type ("none",
 * "manual", or "standardization").
 * @param scaling_type A string specifying the normalization type ("none",
 * "manual", or "standardization").
 * @param modify_x If true, modifies x in-place; otherwise, x remains unchanged
 * (centers/scales are still computed).
 * @param observation_weights Optional nonnegative weights, one per row.
 *
 * @return true if normalization succeeds, false otherwise.
 */
template<typename T>
JitNormalization
normalize(Eigen::MatrixBase<T>& x,
          Eigen::VectorXd& x_centers,
          Eigen::VectorXd& x_scales,
          const std::string& centering_type,
          const std::string& scaling_type,
          const bool modify_x,
          const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
{
  const int p = x.cols();

  computeCenters(x_centers, x, centering_type, observation_weights);
  computeScales(x_scales, x, scaling_type, observation_weights);

  // Handle zero-variance columns by setting their scale to 1.0
  // to avoid division by zero during normalization
  if (scaling_type != "none" && scaling_type != "manual") {
    for (int j = 0; j < p; ++j) {
      if (std::abs(x_scales(j)) == 0.0) {
        x_scales(j) = 1.0;
      }
    }
  }

  bool center = centering_type != "none";
  bool scale = scaling_type != "none";
  bool center_jit = center && !modify_x;
  bool scale_jit = scale && !modify_x;

  JitNormalization jit_normalization;

  if (center_jit && scale_jit) {
    jit_normalization = JitNormalization::Both;
  } else if (center_jit) {
    jit_normalization = JitNormalization::Center;
  } else if (scale_jit) {
    jit_normalization = JitNormalization::Scale;
  } else {
    jit_normalization = JitNormalization::None;
  }

  if (modify_x && (center || scale)) {
    const bool ignored_rows = observation_weights.size() > 0 &&
                              (observation_weights.array() == 0.0).any();
    if (ignored_rows) {
      for (int j = 0; j < p; ++j) {
        for (Eigen::Index i = 0; i < x.rows(); ++i) {
          // Ignored extremes must not overflow and contaminate later products.
          if (observation_weights(i) == 0.0) {
            x(i, j) = 0.0;
            continue;
          }
          if (center) {
            x(i, j) -= x_centers(j);
          }
          if (scale) {
            x(i, j) /= x_scales(j);
          }
        }
      }
    } else {
      for (int j = 0; j < p; ++j) {
        if (center) {
          x.col(j).array() -= x_centers(j);
        }
        if (scale) {
          x.col(j).array() /= x_scales(j);
        }
      }
    }
  }

  return jit_normalization;
}

/**
 * Normalize a sparse matrix by scaling only.
 *
 * To preserve sparsity, centering is not applied. The scaling factors for each
 * column are computed according to the specified normalization type ("none",
 * "manual", or "standardization"). If modify_x is true, the scaling is applied
 * directly to the input matrix.
 *
 * @param x The sparse input matrix to be normalized.
 * @param x_centers A vector that will hold the column centers.
 * For sparse matrices, centering is typically skipped; this
 * parameter is maintained for consistency.
 * @param x_scales  A vector that will hold the column scaling factors. It will
 * be resized to match the number of columns.
 * @param centering_type A string specifying the normalization type ("none",
 * "manual", or "standardization").
 * @param scaling_type A string specifying the normalization type ("none",
 * "manual", or "standardization").
 * @param observation_weights Optional nonnegative weights, one per row.
 *
 * @return true if normalization succeeds, false otherwise.
 */
template<typename T>
JitNormalization
normalize(Eigen::SparseMatrixBase<T>& x,
          Eigen::VectorXd& x_centers,
          Eigen::VectorXd& x_scales,
          const std::string& centering_type,
          const std::string& scaling_type,
          const bool,
          const Eigen::VectorXd& observation_weights = Eigen::VectorXd())
{
  computeCenters(x_centers, x, centering_type, observation_weights);
  computeScales(x_scales, x, scaling_type, observation_weights);

  // Handle zero-variance columns by setting their scale to 1.0
  // to avoid division by zero during normalization
  const int p = x.cols();
  if (scaling_type != "none" && scaling_type != "manual") {
    for (int j = 0; j < p; ++j) {
      if (std::abs(x_scales(j)) == 0.0) {
        x_scales(j) = 1.0;
      }
    }
  }

  bool center = centering_type != "none";
  bool scale = scaling_type != "none";
  bool center_jit = center;
  bool scale_jit = scale;

  JitNormalization jit_normalization;

  if (center_jit && scale_jit) {
    jit_normalization = JitNormalization::Both;
  } else if (center_jit) {
    jit_normalization = JitNormalization::Center;
  } else if (scale_jit) {
    jit_normalization = JitNormalization::Scale;
  } else {
    jit_normalization = JitNormalization::None;
  }

  // TODO: Implement in-place scaling for sparse matrices.
  // if (modify_x && scaling_type != "none") {
  //   for (int j = 0; j < x.cols(); ++j) {
  //     for (Eigen::SparseMatrix<double>::InnerIterator it(x, j); it; ++it) {
  //       it.valueRef() = it.value() / x_scales(j);
  //     }
  //   }
  // }

  return jit_normalization;
}

/**
 * @brief Rescales the coefficients using the given parameters.
 *
 * This function rescales the coefficients by dividing each coefficient by the
 * corresponding scale factor and subtracting the product of the center and
 * the coefficient from the intercept.
 *
 * @param beta0 The intercept coefficient.
 * @param beta The vector of coefficients.
 * @param x_centers The vector of center values.
 * @param x_scales The vector of scale factors.
 * @param intercept Should an intercept be fit?
 * @return A tuple containing the rescaled intercept and coefficients.
 *
 * @note The input vectors `beta`, `x_centers`, and `x_scales` must have the
 * same size.
 * @note The output vector `beta` will be modified in-place.
 *
 * @see SlopeParameters
 */
std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
rescaleCoefficients(const Eigen::VectorXd& beta0,
                    const Eigen::SparseMatrix<double>& beta,
                    const Eigen::VectorXd& x_centers,
                    const Eigen::VectorXd& x_scales,
                    const bool intercept);

} // namespace slope
