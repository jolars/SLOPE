/**
 * @internal
 * @file
 * @brief Mathematical support functions for the slope package.
 */

#pragma once

#include "jit_normalization.h"
#include "threads.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <numeric>
#include <vector>

namespace slope {

/**
 * @brief Returns the sign of a given value.
 *
 * This function determines the sign of the input value by comparing it to zero.
 * It returns -1 if the value is negative, 0 if the value is zero, and 1 if the
 * value is positive.
 *
 * @tparam T The type of the input value.
 * @param val The value for which the sign needs to be determined.
 * @return -1 if the value is negative, 0 if the value is zero, and 1 if the
 * value is positive.
 */
template<typename T>
int
sign(T val)
{
  return (T(0) < val) - (val < T(0));
}

/**
 * Calculates the cumulative sum of the elements in the input array.
 *
 * @tparam T The type of the input array.
 * @param x The input array.
 * @return An Eigen::ArrayXd containing the cumulative sum of the elements in
 * the input array.
 */
template<typename T>
Eigen::ArrayXd
cumSum(const T& x)
{
  std::vector<double> cum_sum(x.size());
  std::partial_sum(
    x.data(), x.data() + x.size(), cum_sum.begin(), std::plus<double>());

  Eigen::Map<Eigen::ArrayXd> out(cum_sum.data(), cum_sum.size());

  return out;
}

/**
 * Calculates the sigmoid function for the given input.
 *
 * The sigmoid function is defined as 1 / (1 + exp(-x)).
 *
 * @tparam T The type of the input.
 * @param x The input value.
 * @return The result of the sigmoid function.
 */
template<typename T>
T
sigmoid(const T& x)
{
  return 1.0 / (1.0 + std::exp(-x));
}

/**
 * The logit
 *
 * The logit function is defined as \f$\log(\frac{x}{1 - x})\f$.
 *
 * @tparam T The type of the input.
 * @param x The input value.
 * @return The result of the logit function.
 */
template<typename T>
T
logit(const T& x)
{
  assert(x > 0 && x < 1 && "Input must be in (0, 1)");

  return std::log(x) - std::log1p(-x);
}

/**
 * Returns the value of x clamped between the specified lower and upper bounds.
 *
 * @tparam T the type of the values being clamped
 * @param x the value to be clamped
 * @param lo the lower bound
 * @param hi the upper bound
 * @return the clamped value of x
 */
template<typename T>
T
clamp(const T& x, const T& lo, const T& hi)
{
  return x < lo ? lo : x > hi ? hi : x;
}

/**
 * LogSumExp
 *
 * @param a A matrix
 * @return \f$\log(\sum_i \exp(a_i))\f$
 */
Eigen::VectorXd
logSumExp(const Eigen::MatrixXd& a);

/**
 * Softmax
 *
 * Computes the softmax function for the given input matrix.
 *
 * @param x A matrix
 * @return \f$\exp(a) / \sum_i \exp(a_i)\f$
 */
Eigen::MatrixXd
softmax(const Eigen::MatrixXd& x);

/**
 * Computes the gradient of the loss with respect to \f(\beta\f).
 *
 * @tparam T The type of the input matrix.
 * @param x The input matrix.
 * @param active_set Indicies for active set
 * @param beta0 Intercept
 * @param beta Coefficients
 * @param x_centers The vector of center values for each column of x.
 * @param x_scales The vector of scale values for each column of x.
 * @param jit_normalization Type of JIT normalization.
 * @param intercept Whether to fit an intercept.
 * @return The computed gradient vector.
 */
template<typename T>
Eigen::MatrixXd
linearPredictor(const T& x,
                const std::vector<int>& active_set,
                const Eigen::VectorXd& beta0,
                const Eigen::VectorXd& beta,
                const Eigen::VectorXd& x_centers,
                const Eigen::VectorXd& x_scales,
                const JitNormalization jit_normalization,
                const bool intercept)
{
  int n = x.rows();
  int p = x.cols();
  int m = beta0.size();

  Eigen::MatrixXd eta = Eigen::MatrixXd::Zero(n, m);

  for (const auto& ind : active_set) {
    auto [k, j] = std::div(ind, p);

    switch (jit_normalization) {
      case JitNormalization::Both:
        eta.col(k) += x.col(j) * beta(ind) / x_scales(j);
        eta.col(k).array() -= beta(ind) * x_centers(j) / x_scales(j);
        break;

      case JitNormalization::Center:
        eta.col(k) += x.col(j) * beta(ind);
        eta.col(k).array() -= beta(ind) * x_centers(j);
        break;

      case JitNormalization::Scale:
        eta.col(k) += x.col(j) * beta(ind) / x_scales(j);
        break;

      case JitNormalization::None:
        eta.col(k) += x.col(j) * beta(ind);
        break;
    }
  }

  if (intercept) {
    eta.rowwise() += beta0.transpose();
  }

  return eta;
}

/**
 * Computes the gradient of the loss with respect to \f(\beta\f).
 *
 * @tparam T The type of the input matrix.
 * @param gradient The gradient vector.
 * @param x The input matrix.
 * @param residual The residual matrix.
 * @param active_set The indices for the active set.
 * @param x_centers The vector of center values for each column of x.
 * @param x_scales The vector of scale values for each column of x.
 * @param w Working weights
 * @param jit_normalization Type of JIT normalization
 * just-in-time.
 */
template<typename T>
void
updateGradient(Eigen::VectorXd& gradient,
               const T& x,
               const Eigen::MatrixXd& residual,
               const std::vector<int>& active_set,
               const Eigen::VectorXd& x_centers,
               const Eigen::VectorXd& x_scales,
               const Eigen::VectorXd& w,
               const JitNormalization jit_normalization)
{
  const int n = x.rows();
  const int p = x.cols();
  const int m = residual.cols();

  assert(gradient.size() == p * m &&
         "Gradient matrix has incorrect dimensions");

  Eigen::MatrixXd weighted_residual(n, m);
  Eigen::ArrayXd wr_sums(m);

#ifdef _OPENMP
#pragma omp parallel for num_threads(Threads::get())
#endif
  for (int k = 0; k < m; ++k) {
    weighted_residual.col(k) = residual.col(k).cwiseProduct(w);
    wr_sums(k) = weighted_residual.col(k).sum();
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(Threads::get())
#endif
  for (size_t i = 0; i < active_set.size(); ++i) {
    int ind = active_set[i];
    auto [k, j] = std::div(ind, p);

    switch (jit_normalization) {
      case JitNormalization::Both:
        gradient(ind) =
          (x.col(j).dot(weighted_residual.col(k)) - x_centers(j) * wr_sums(k)) /
          (x_scales(j) * n);
        break;
      case JitNormalization::Center:
        gradient(ind) =
          (x.col(j).dot(weighted_residual.col(k)) - x_centers(j) * wr_sums(k)) /
          n;
        break;
      case JitNormalization::Scale:
        gradient(ind) =
          x.col(j).dot(weighted_residual.col(k)) / (x_scales(j) * n);
        break;
      case JitNormalization::None:
        gradient(ind) = x.col(j).dot(weighted_residual.col(k)) / n;
        break;
    }
  }
}

/**
 * Computes the gradient of the loss with respect to \f(\beta\f).
 *
 * @tparam T The type of the input matrix.
 * @param gradient The residual vector.
 * @param x The input matrix.
 * @param offset Gradient offset
 * @param active_set Indices for the active_set
 * @param x_centers The vector of center values for each column of x.
 * @param x_scales The vector of scale values for each column of x.
 * @param jit_normalization Type of JIT normalization
 * just-in-time.
 */
template<typename T>
void
offsetGradient(Eigen::VectorXd& gradient,
               const T& x,
               const Eigen::VectorXd& offset,
               const std::vector<int>& active_set,
               const Eigen::VectorXd& x_centers,
               const Eigen::VectorXd& x_scales,
               const JitNormalization jit_normalization)
{
  const int n = x.rows();
  const int p = x.cols();

  for (size_t i = 0; i < active_set.size(); ++i) {
    int ind = active_set[i];
    auto [k, j] = std::div(ind, p);

    switch (jit_normalization) {
      case JitNormalization::Both:
        gradient(ind) -=
          offset(k) * (x.col(j).sum() / n - x_centers(j)) / x_scales(j);
        break;
      case JitNormalization::Center:
        gradient(ind) -= offset(k) * (x.col(j).sum() / n - x_centers(j));
        break;
      case JitNormalization::Scale:
        gradient(ind) -= offset(k) * x.col(j).sum() / (n * x_scales(j));
        break;
      case JitNormalization::None:
        gradient(ind) -= offset(k) * x.col(j).sum() / n;
        break;
    }
  }
}

/**
 * @brief Computes the union of two sorted integer vectors
 *
 * @param a First sorted vector of integers
 * @param b Second sorted vector of integers
 * @return std::vector<int> Vector containing all elements that appear in either
 * a or b, without duplicates and in sorted order
 */
std::vector<int>
setUnion(const std::vector<int>& a, const std::vector<int>& b);

/**
 * @brief Computes the set difference of two sorted integer vectors
 *
 * @param a First sorted vector of integers (set to subtract from)
 * @param b Second sorted vector of integers (set to subtract)
 * @return std::vector<int> Vector containing elements in a that do not appear
 * in b, maintaining sorted order
 *
 * Returns A \ B = {x ∈ A | x ∉ B}
 */
std::vector<int>
setDiff(const std::vector<int>& a, const std::vector<int>& b);

/**
 * @brief Returns the index of the maximum element in a container
 *
 * @tparam T Container type that supports iterators and std::max_element
 * @param x Container whose maximum element's index is to be found
 * @return int Zero-based index position of the maximum element
 *
 * Uses std::max_element to find the iterator to the maximum element,
 * then converts to index position using std::distance.
 * For containers with multiple maximum elements, returns the first occurrence.
 */
template<typename T>
int
whichMax(const T& x)
{
  return std::distance(x.begin(), std::max_element(x.begin(), x.end()));
}

/**
 * @brief Returns the index of the minimum element in a container
 *
 * @tparam T Container type that supports iterators and std::min_element
 * @param x Container whose maximum element's index is to be found
 * @return int Zero-based index position of the maximum element
 *
 * Uses std::max_element to find the iterator to the maximum element,
 * then converts to index position using std::distance.
 * For containers with multiple maximum elements, returns the first occurrence.
 */
template<typename T>
int
whichMin(const T& x)
{
  return std::distance(x.begin(), std::min_element(x.begin(), x.end()));
}

/**
 * @brief Returns the index of the minimum element in a container
 *
 * @tparam T Container type that supports iterators and std::min_element
 * @param x Container whose maximum element's index is to be found
 * @param comp Comparator that returns true if the first argument is worse
 * than the second argument
 * @return int Zero-based index position of the maximum element
 *
 * Uses std::max_element to find the iterator to the maximum element,
 * then converts to index position using std::distance.
 * For containers with multiple maximum elements, returns the first occurrence.
 */
template<typename T, typename Comparator>
int
whichBest(const T& x, const Comparator& comp)
{
  return std::distance(x.begin(), std::max_element(x.begin(), x.end(), comp));
}

/**
 * @brief Creates an array of n numbers in geometric progression from start to
 * end.
 *
 * @param start The starting value of the sequence
 * @param end The final value of the sequence
 * @param n The number of points to generate
 * @return Eigen::ArrayXd Array containing n points spaced geometrically between
 * start and end
 *
 * @note Similar to numpy.geomspace, generates points that are evenly spaced on
 * a log scale
 * @throws std::invalid_argument If start or end is zero, or if n < 1
 */
Eigen::ArrayXd
geomSpace(const double start, const double end, const int n);

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

} // namespace slope
