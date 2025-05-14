/**
 * @file
 * @brief Various utility functions
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <algorithm>
#include <numeric>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

namespace slope {

/**
 * Sorts the elements in a container in ascending or descending order.
 *
 * This function sorts the elements in the container `v` in either ascending or
 * descending order, depending on the value of the `descending` parameter.
 *
 * @tparam T The type of the container.
 * @param v The container to be sorted.
 * @param descending Flag indicating whether to sort in descending order
 * (default is ascending order).
 *
 * @note The elements in the container must be comparable using the `<` and `>`
 * operators.
 *
 * @see std::sort
 */
template<typename T>
void
sort(T& v, const bool descending = false)
{
  if (descending) {
    std::sort(v.data(), v.data() + v.size(), std::greater<double>());
  } else {
    std::sort(v.data(), v.data() + v.size(), std::less<double>());
  }
}

/**
 * Returns indices of true values in a boolean container.
 *
 * @tparam T Container type supporting size() and operator[] (e.g.,
 * std::vector<bool>, std::array<bool>)
 * @param x Input container with boolean-convertible values
 * @return std::vector<int> containing indices where x[i] evaluates to true
 *
 * Example:
 *   std::vector<bool> v = {true, false, true, false, true};
 *   auto indices = which(v); // returns {0, 2, 4}
 */
template<typename T>
std::vector<int>
which(const T& x)
{
  std::vector<int> out;
  for (int i = 0; i < x.size(); i++) {
    if (x[i]) {
      out.emplace_back(i);
    }
  }

  return out;
}

/**
 * Sorts the elements of a vector and returns the indices of the sorted
 * elements.
 *
 * This function sorts the elements of a vector in ascending or descending order
 * and returns the indices of the sorted elements. The sorting is done using the
 * std::sort function from the C++ standard library.
 *
 * @tparam T The type of the vector elements.
 * @param v The vector to be sorted.
 * @param descending Flag indicating whether to sort in descending order.
 * Default is false.
 * @return A vector of indices representing the sorted order of the elements in
 * the input vector.
 */
template<typename T>
std::vector<int>
sortIndex(T& v, const bool descending = false)
{
  using namespace std;

  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  if (descending) {
    sort(idx.begin(), idx.end(), [&v](int i, int j) { return v[i] > v[j]; });
  } else {
    sort(idx.begin(), idx.end(), [&v](int i, int j) { return v[i] < v[j]; });
  }

  return idx;
}

/**
 * Permutes the elements of a container according to the given indices.
 *
 * This function takes a container of values and permutes its elements according
 * to the given indices. The indices specify the new order of the elements in
 * the container.
 *
 * @tparam T The type of the container.
 * @param values The container of values to be permuted.
 * @param ind The vector of indices specifying the new order of the elements.
 */
template<typename T>
void
permute(T& values, const std::vector<int>& ind)
{
  /**
   * @brief The container to store the permuted values.
   */
  T out(values.size());

  /**
   * @brief Permute the values according to the given indices.
   */
  for (int i = 0; i < values.size(); ++i)
    out[i] = std::move(values[ind[i]]);

  /**
   * @brief Assign the permuted values back to the original container.
   */
  values = std::move(out);
}

/**
 * Inverse permutes the elements of a container based on the given
 * indices.
 *
 * This function takes a container of values and a vector of indices and
 * rearranges the elements of the container according to the indices. The
 * resulting container will have the elements in the order specified by the
 * indices.
 *
 * @tparam T The type of the container.
 * @param values The container of values to be permuted.
 * @param ind The vector of indices specifying the new order of the elements.
 */
template<typename T>
void
inversePermute(T& values, const std::vector<int>& ind)
{
  T out(values.size()); /**< The resulting container after permutation. */

  for (int i = 0; i < values.size(); ++i)
    out[ind[i]] = std::move(values[i]);

  values = std::move(out);
}

/**
 * Moves a range of elements within a vector.
 *
 * This function moves a range of elements within a vector from one position to
 * another. The elements are moved in a way that the order is preserved.
 *
 * @tparam T The type of elements in the vector.
 * @param v The vector containing the elements.
 * @param from The starting index of the range to be moved.
 * @param to The ending index of the range to be moved.
 * @param size The size of the range to be moved.
 */
template<typename T>
void
move_elements(std::vector<T>& v, const int from, const int to, const int size)
{
  if (from > to) {
    std::rotate(v.begin() + to, v.begin() + from, v.begin() + from + size);
  } else {
    std::rotate(v.begin() + from, v.begin() + from + size, v.begin() + to + 1);
  }
}

/**
 * @brief Validates if a given value exists in a set of valid options
 *
 * @details Throws an informative error if the value is not found in the valid
 * options, listing all valid possibilities in the error message.
 *
 * @param value The value to validate
 * @param valid_options Set of valid options
 * @param parameter_name Name of the parameter being validated (used in error
 * message).
 * @throws std::invalid_argument If value is not in valid_options
 */
void
validateOption(const std::string& value,
               const std::set<std::string>& valid_options,
               const std::string& parameter_name);

/**
 * @brief Extract a subset of rows from an Eigen matrix
 *
 * @param x The input matrix to extract rows from
 * @param indices A vector of row indices to extract
 * @return Eigen::MatrixXd A new matrix containing only the specified rows
 *
 * This function creates a new matrix containing only the rows specified in the
 * indices vector, preserving their order. The number of columns remains the
 * same.
 */
Eigen::MatrixXd
subset(const Eigen::MatrixXd& x, const std::vector<int>& indices);

/**
 * @brief Extract a subset of rows from a sparse Eigen matrix
 *
 * @param x The input sparse matrix to extract rows from
 * @param indices A vector of row indices to extract
 * @return Eigen::SparseMatrix<double> A new sparse matrix containing only the
 * specified rows
 *
 * This function creates a new sparse matrix containing only the rows specified
 * in the indices vector, preserving their order. The number of columns remains
 * the same. The sparsity structure is maintained in the extracted rows.
 */
Eigen::SparseMatrix<double>
subset(const Eigen::SparseMatrix<double>& x, const std::vector<int>& indices);

/**
 * @brief Extract specified columns from a dense matrix
 *
 * Creates a new matrix containing only the columns specified in the indices
 * vector.
 *
 * @tparam MatrixType The type of matrix (supports both Eigen::MatrixXd and
 * Eigen::Map<>)
 * @param x Input matrix
 * @param indices Vector of column indices to extract (0-based)
 * @return Eigen::MatrixXd New matrix with only the selected columns
 */
template<typename MatrixType>
Eigen::MatrixXd
subsetCols(const MatrixType& x, const std::vector<int>& indices)
{
  if (indices.empty()) {
    return Eigen::MatrixXd(x.rows(), 0);
  }

  Eigen::MatrixXd result(x.rows(), indices.size());
  for (size_t i = 0; i < indices.size(); ++i) {
    result.col(i) = x.col(indices[i]);
  }
  return result;
}

/**
 * @brief Create a set of unique values from an Eigen matrix
 *
 * @param x The input matrix to extract unique values from
 * @return std::unordered_set<double> A set containing all unique values found
 * in the matrix
 *
 * This function iterates through all elements of the input matrix and collects
 * all unique values into an unordered set. The order of elements in the
 * returned set is not guaranteed.
 */
inline std::unordered_set<double>
unique(const Eigen::MatrixXd& x)
{
  std::unordered_set<double> unique;
  for (Eigen::Index j = 0; j < x.cols(); j++) {
    for (Eigen::Index i = 0; i < x.rows(); i++) {
      unique.insert(x(i, j));
    }
  }

  return unique;
}

} // namespace slope
