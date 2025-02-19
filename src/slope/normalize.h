/** @file
 * @brief Functions to normalize the design matrix and rescale coefficients
 * in case the design was normalized
 */

#pragma once

#include <Eigen/SparseCore>

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

/**
 * @brief Computes the L1 (Manhattan) norms for each column of a matrix
 *
 * @tparam T Matrix type (must support cols(), col() operations, compatible with
 * Eigen)
 * @param x Input matrix whose column L1 norms are to be computed
 * @return Eigen::VectorXd Vector containing L1 norms of each column
 *
 * For each column j in matrix x, computes the L1 norm:
 * \f[ \|x_j\|_1 = \sum_{i} |x_{ij}| \f]
 * where x_{ij} represents the i-th element of the j-th column.
 */
template<typename T>
Eigen::VectorXd
l1Norms(const T& x)
{
  const int p = x.cols();

  Eigen::VectorXd out(p);

  for (int j = 0; j < p; ++j) {
    out(j) = x.col(j).cwiseAbs().sum();
  }

  return out;
}

/**
 * @brief Computes the L2 (Euclidean) norms for each column of a sparse matrix
 *
 * @param x Input sparse matrix whose column L2 norms are to be computed
 * @return Eigen::VectorXd Vector containing L2 norms of each column
 *
 * For each column j in matrix x, computes the L2 norm:
 * \f[ \|x_j\|_2 = \sqrt{\sum_{i} x_{ij}^2} \f]
 * where x_{ij} represents the i-th element of the j-th column.
 */
Eigen::VectorXd
l2Norms(const Eigen::SparseMatrix<double>& x);

/**
 * @brief Computes the L2 (Euclidean) norms for each column of a dense matrix
 *
 * @param x Input dense matrix whose column L2 norms are to be computed
 * @return Eigen::VectorXd Vector containing L2 norms of each column
 *
 * For each column j in matrix x, computes the L2 norm:
 * \f[ \|x_j\|_2 = \sqrt{\sum_{i} x_{ij}^2} \f]
 * where x_{ij} represents the i-th element of the j-th column.
 */
Eigen::VectorXd
l2Norms(const Eigen::MatrixXd& x);

/**
 * @brief Computes the maximum absolute value for each column of a matrix
 *
 * @param x Input matrix to compute column-wise maximum absolute values
 * @return Eigen::VectorXd Vector containing maximum absolute values for each
 * column
 *
 * For each column j in matrix x, computes:
 * \f[ \max_{i} |x_{ij}| \f]
 * where x_{ij} represents the i-th element of the j-th column.
 */
Eigen::VectorXd
maxAbs(const Eigen::SparseMatrix<double>& x);

/**
 * @brief Computes the maximum absolute value for each column of a dense matrix
 *
 * @param x Input dense matrix whose column-wise maximum absolute values are to
 * be computed
 * @return Eigen::VectorXd Vector containing the maximum absolute value of each
 * column
 *
 * For each column j in matrix x, computes:
 * \f[ \max_{i} |x_{ij}| \f]
 * where x_{ij} represents the i-th element of the j-th column.
 */
Eigen::VectorXd
maxAbs(const Eigen::MatrixXd& x);

/**
 * @brief Computes the arithmetic mean for each column of a sparse matrix
 *
 * @param x Input sparse matrix whose column means are to be computed
 * @return Eigen::VectorXd Vector containing the arithmetic mean of each column
 *
 * For each column j in matrix x, computes the mean:
 * \f[ \bar{x}_j = \frac{1}{n}\sum_{i} x_{ij} \f]
 * where n is the number of rows in the matrix and x_{ij} represents
 * the i-th element of the j-th column.
 */
Eigen::VectorXd
means(const Eigen::SparseMatrix<double>& x);

/**
 * @brief Computes the arithmetic mean for each column of a dense matrix
 *
 * @param x Input dense matrix whose column means are to be computed
 * @return Eigen::VectorXd Vector containing the arithmetic mean of each column
 *
 * For each column j in matrix x, computes the mean:
 * \f[ \bar{x}_j = \frac{1}{n}\sum_{i} x_{ij} \f]
 * where n is the number of rows in the matrix and x_{ij} represents
 * the i-th element of the j-th column.
 */
Eigen::VectorXd
means(const Eigen::MatrixXd& x);

/**
 * @brief Computes the standard deviation for each column of a matrix
 *
 * @tparam T Matrix type (must support sparse matrix operations)
 * @param x Input matrix whose column standard deviations are to be computed
 * @return Eigen::VectorXd Vector containing the standard deviation of each
 * column
 *
 * For each column j in matrix x, computes the standard deviation:
 * \f[ \sigma_j = \sqrt{\frac{1}{n}\sum_{i} (x_{ij} - \bar{x}_j)^2} \f]
 * where n is the number of rows, x_{ij} represents the i-th element of the j-th
 * column, and \f$\bar{x}_j\f$ is the mean of column j.
 *
 * This function uses Welford's algorithm.
 */
template<typename T>
Eigen::VectorXd
stdDevs(const T& x)
{
  const int n = x.rows();
  const int p = x.cols();

  Eigen::VectorXd x_means = means(x);
  Eigen::VectorXd out(p);

  for (int j = 0; j < p; ++j) {
    double m2 = 0.0;
    const double mean = x_means(j);

    // Process non-zero elements
    for (typename T::InnerIterator it(x, j); it; ++it) {
      double delta = it.value() - mean;
      m2 += delta * delta;
    }

    // Account for zeros
    int nz_count = x.col(j).nonZeros();
    if (nz_count < n) {
      m2 += (n - nz_count) * mean * mean;
    }

    out(j) = std::sqrt(m2 / n);
  }

  return out;
}

/**
 * @brief Computes the range (max - min) for each column of a matrix
 *
 * @param x Input matrix whose column ranges are to be computed
 * @return Eigen::VectorXd Vector containing the range of each column
 *
 * For each column j in matrix x, computes:
 * \f[ range_j = \max_{i}(x_{ij}) - \min_{i}(x_{ij}) \f]
 * where x_{ij} represents the i-th element of the j-th column.
 */
Eigen::VectorXd
ranges(const Eigen::SparseMatrix<double>& x);

/**
 * @brief Computes the range (max - min) for each column of a dense matrix
 *
 * @param x Input dense matrix whose column ranges are to be computed
 * @return Eigen::VectorXd Vector containing the range of each column
 *
 * For each column j in matrix x, computes:
 * \f[ range_j = \max_{i}(x_{ij}) - \min_{i}(x_{ij}) \f]
 * where x_{ij} represents the i-th element of the j-th column.
 *
 * Uses Eigen's efficient colwise() operations to compute
 * maximum and minimum values for each column simultaneously.
 */
Eigen::VectorXd
ranges(const Eigen::MatrixXd& x);

/**
 * @brief Computes the minimum value for each column of a sparse matrix
 *
 * @param x Input sparse matrix whose column minimums are to be computed
 * @return Eigen::VectorXd Vector containing the minimum value of each column
 *
 * For each column j in matrix x, computes:
 * \f[ \min_{i}(x_{ij}) \f]
 * where x_{ij} represents the i-th element of the j-th column.
 */
Eigen::VectorXd
mins(const Eigen::SparseMatrix<double>& x);

/**
 * @brief Computes the minimum value for each column of a dense matrix
 *
 * @param x Input dense matrix whose column minimums are to be computed
 * @return Eigen::VectorXd Vector containing the minimum value of each column
 *
 * For each column j in matrix x, computes:
 * \f[ \min_{i}(x_{ij}) \f]
 * where x_{ij} represents the i-th element of the j-th column.
 *
 * Uses Eigen's built-in column-wise operations for efficient computation
 * on dense matrices.
 */
Eigen::VectorXd
mins(const Eigen::MatrixXd& x);

/**
 * Compute centers.
 *
 * There are two supported centering types:
 * - "none": Do not compute centers and scales.
 * - "mean": Use arithmetic means.
 *
 * @tparam T The type of the input matrix.
 * @param x_centers A vector where the computed or provided centers will be
 * stored.
 * @param x The input matrix.
 * @param type A string specifying the normalization type ("none", "manual", or
 * "standardization").
 *
 * @throws std::invalid_argument if the provided manual centers or scales have
 * invalid dimensions or contain non-finite values.
 */
template<typename T>
void
computeCenters(Eigen::VectorXd& x_centers, const T& x, const std::string& type)
{
  int p = x.cols();

  if (type == "manual") {
    if (x_centers.size() != p) {
      throw std::invalid_argument("Invalid dimensions in centers");
    }

    if (!x_centers.allFinite()) {
      throw std::invalid_argument("Centers must be finite");
    }

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
 * There are two supported scaling types:
 * - "none": Do not compute centers and scales.
 * - "sd": Compute centers and scales using Welford’s algorithm.
 *
 * @tparam T The type of the input matrix.
 * @param x The input matrix.
 * @param x_scales A vector where the computed or provided scales will be
 * stored.
 * @param type A string specifying the normalization type ("none", "manual", or
 * "standardization").
 *
 * @throws std::invalid_argument if the provided manual centers or scales have
 * invalid dimensions or contain non-finite values.
 */
template<typename T>
void
computeScales(Eigen::VectorXd& x_scales, const T& x, const std::string& type)
{
  int p = x.cols();

  if (type == "manual") {
    if (x_scales.size() != p) {
      throw std::invalid_argument("Invalid dimensions in scales");
    }
    if (!x_scales.allFinite()) {
      throw std::invalid_argument("Scales must be finite");
    }
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
 *
 * @return true if normalization succeeds, false otherwise.
 */
JitNormalization
normalize(Eigen::MatrixXd& x,
          Eigen::VectorXd& x_centers,
          Eigen::VectorXd& x_scales,
          const std::string& centering_type,
          const std::string& scaling_type,
          const bool modify_x);

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
 * @param modify_x If true, performs in-place scaling on x; otherwise, leaves x
 * unchanged.
 *
 * @return true if normalization succeeds, false otherwise.
 */
JitNormalization
normalize(Eigen::SparseMatrix<double>& x,
          Eigen::VectorXd& x_centers,
          Eigen::VectorXd& x_scales,
          const std::string& centering_type,
          const std::string& scaling_type,
          const bool modify_x);

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
rescaleCoefficients(Eigen::VectorXd beta0,
                    Eigen::MatrixXd beta,
                    const Eigen::VectorXd& x_centers,
                    const Eigen::VectorXd& x_scales,
                    const bool intercept);

} // namespace slope
