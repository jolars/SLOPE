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
  int n_clusters() const;

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
  std::vector<double> coeffs() const;

  /**
   * @brief Returns a vector containing the indices of all clusters.
   * @return A vector containing the indices.
   */
  std::vector<int> indices() const;

  /**
   * @brief Returns a vector containing the pointers of all clusters.
   * @return A vector containing the pointers.
   */
  std::vector<int> pointers() const;

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

  /**
   * @brief Returns the cluster pattern as a sparse matrix
   * @return A sparse matrix, where the indices of column j
   *   are the indices of the features in cluster j.
   *   The matrix is of size \f$p \times (m - 1)\f$, where \f$p\f$ is the number
   *   of features and \f$m\f$ is the number of clusters. The features
   *   of the zero cluster are not included.
   */
  Eigen::SparseMatrix<int> patternMatrix() const;

private:
  std::vector<double> c;  /**< The coefficients of the clusters. */
  std::vector<int> c_ind; /**< The indices of the clusters. */
  std::vector<int> signs; /**< The signs of the coefficients. */
  std::vector<int>
    c_ptr; /**< Pointers to the start of each of the clusters' indices. */
  int p;   /**< The number of features. */

  /**
   * @brief Mutable vector to store zero indices when needed
   * This is mutable to allow const methods to modify it
   */
  mutable std::vector<int> zero_indices;

  /**
   * @brief Flag to track if zero_indices are up to date
   */
  mutable bool zero_indices_valid = false;

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

  /**
   * @brief Helper function to compute zero cluster indices
   * @return Reference to vector containing zero cluster indices
   */
  std::vector<int>& getZeroIndices() const;

  /**
   * @brief Checks if this represents an all-zeros vector
   * @return true if the cluster contains only zeros
   */
  bool hasAllZeros() const;
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
