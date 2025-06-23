/**
 * @file
 * @brief The declaration of the Clusters class
 */

#pragma once

#include <Eigen/SparseCore>
#include <vector>

namespace slope {

/**
 * @class Clusters
 * @brief Representation of the clusters in SLOPE
 */
class Clusters
{
public:
  /**
   * @brief Constructs an Clusters object.
   */
  Clusters() = default;

  /**
   * @brief Constructs a Clusters object with the given beta vector.
   * @param beta The beta vector.
   */
  Clusters(const Eigen::VectorXd& beta);

  /**
   * @brief Returns the number of clusters.
   * @return The size of the cluster.
   */
  std::size_t size() { return c.size(); }

  /**
   * @brief Returns an iterator pointing to the beginning of the cluster with
   * the given index.
   * @param i The index of the cluster.
   * @return An iterator pointing to the beginning of the cluster.
   */
  std::vector<int>::iterator begin(const int i);

  /**
   * @brief Returns an iterator pointing to the end of the cluster with the
   * given index.
   * @param i The index of the cluster.
   * @return An iterator pointing to the end of the cluster.
   */
  std::vector<int>::iterator end(const int i);

  /**
   * @brief Returns a constant iterator pointing to the beginning of the cluster
   * with the given index.
   * @param i The index of the cluster.
   * @return A constant iterator pointing to the beginning of the cluster.
   */
  std::vector<int>::const_iterator cbegin(const int i) const;

  /**
   * @brief Returns a constant iterator pointing to the end of the cluster with
   * the given index.
   * @param i The index of the cluster.
   * @return A constant iterator pointing to the end of the cluster.
   */
  std::vector<int>::const_iterator cend(const int i) const;

  /**
   * @brief Returns the size of the cluster with the given index.
   * @param i The index of the cluster.
   * @return The size of the cluster.
   */
  int cluster_size(const int i) const;

  /**
   * @brief Returns the pointer of the cluster with the given index.
   * @param i The index of the cluster.
   * @return The pointer of the cluster.
   */
  int pointer(const int i) const;

  /**
   * @brief Returns the number of clusters.
   * @return The number of clusters.
   */
  int size() const;

  /**
   * @brief Returns the coefficient of the cluster with the given index.
   * @param i The index of the cluster.
   * @return The coefficient of the cluster.
   */
  double coeff(const int i) const;

  /**
   * @brief Sets the coefficient of the cluster with the given index.
   * @param i The index of the cluster.
   * @param x The new coefficient value.
   */
  void setCoeff(const int i, const double x);

  /**
   * @brief Returns a vector containing the coefficients of all clusters.
   * @return A vector containing the coefficients.
   */
  const std::vector<double>& coeffs() const;

  /**
   * @brief Returns a vector containing the indices of all clusters.
   * @return A vector containing the indices.
   */
  const std::vector<int>& indices() const;

  /**
   * @brief Returns a vector containing the pointers of all clusters.
   * @return A vector containing the pointers.
   */
  const std::vector<int>& pointers() const;

  /**
   * @brief Updates the cluster structure when an index is changed.
   * @param old_index The old index.
   * @param new_index The new index.
   * @param c_new The new coefficient value.
   */
  void update(const int old_index, const int new_index, const double c_new);

  /**
   * @brief Updates the cluster structure with the given beta vector.
   * @param beta The beta vector.
   */
  void update(const Eigen::VectorXd& beta);

  /**
   * @brief Returns the clusters as a vector of vectors.
   * @return The clusters as a vector of vectors.
   */
  std::vector<std::vector<int>> getClusters() const;

private:
  std::vector<double> c;  /**< The coefficients of the clusters. */
  std::vector<int> c_ind; /**< The indices of the clusters. */
  std::vector<int>
    c_ptr; /**< Pointers to the start of each of the clusters' indices. */
  int p;   /**< The number of features. */

  /**
   * @brief Reorders the cluster structure when an index is changed.
   * @param old_index The old index.
   * @param new_index The new index.
   */
  void reorder(const int old_index, const int new_index);

  /**
   * @brief Merges two clusters into one.
   * @param old_index The index of the cluster to be merged.
   * @param new_index The index of the cluster to merge into.
   */
  void merge(const int old_index, const int new_index);
};

/**
 * @brief Returns the cluster pattern as a sparse matrix
 * @param beta The beta vector used to create the clusters
 * @return A sparse matrix, where the indices of column j
 *   are the indices of the features in cluster j.
 *   The matrix is of size \f$p \times (m - 1)\f$, where \f$p\f$ is the number
 *   of features and \f$m\f$ is the number of clusters. The features
 *   of the zero cluster are not included.
 */
Eigen::SparseMatrix<int>
patternMatrix(const Eigen::VectorXd& beta);

} // namespace slope
