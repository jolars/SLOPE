/**
 * @file
 * @brief Cross-validation fold management for SLOPE models
 *
 * This file provides functionality for creating and managing data partitions
 * for k-fold cross-validation. It supports both automatic fold generation
 * with random partitioning and user-defined custom folds.
 */
#pragma once

#include "utils.h"
#include <Eigen/Core>
#include <cstdint>
#include <vector>

namespace slope {

/**
 * @brief Manages data partitioning for cross-validation
 *
 * This class handles the creation and access of training/test splits
 * for k-fold cross-validation. It can generate random folds or use
 * user-provided custom fold definitions.
 */
class Folds
{
public:
  /**
   * @brief Constructor for generating random folds
   *
   * @param n_samples Total number of samples in the dataset
   * @param n_folds Number of folds to create
   * @param seed Random seed for reproducibility
   *
   * Creates a k-fold partition by randomly assigning samples to folds.
   */
  Folds(int n_samples, int n_folds, uint64_t seed = 42)
    : folds(createFolds(n_samples, n_folds, seed))
    , n_folds(n_folds)
  {
  }

  /**
   * @brief Constructor for user-provided folds
   *
   * @param folds Vector of vectors containing indices for each fold
   *
   * Uses the provided fold assignments instead of creating random folds.
   */
  explicit Folds(std::vector<std::vector<int>> folds)
    : folds(std::move(folds))
    , n_folds(folds.size())
  {
  }

  /**
   * @brief Get test indices for a specific fold
   *
   * @param fold_idx Index of the fold to use as test set
   * @return const std::vector<int>& Reference to the vector of test indices
   */
  const std::vector<int>& getTestIndices(size_t fold_idx) const;

  /**
   * @brief Get training indices for a specific fold
   *
   * @param fold_idx Index of the fold to exclude from training
   * @return std::vector<int> Vector containing all indices except those in the
   * specified fold
   */
  std::vector<int> getTrainingIndices(size_t fold_idx) const;

  /**
   * @brief Split data into training and test sets for a specific fold
   *
   * @tparam MatrixType Type of design matrix (supports both dense and sparse
   * matrices)
   * @param x Input feature matrix
   * @param y Response matrix
   * @param fold_idx Index of the fold to use as test set
   * @return std::tuple containing (x_train, y_train, x_test, y_test)
   */
  template<typename MatrixType>
  std::tuple<MatrixType, Eigen::MatrixXd, MatrixType, Eigen::MatrixXd>
  split(MatrixType& x, const Eigen::MatrixXd& y, size_t fold_idx) const
  {
    auto test_idx = getTestIndices(fold_idx);
    auto train_idx = getTrainingIndices(fold_idx);

    MatrixType x_test = subset(x, test_idx);
    Eigen::MatrixXd y_test = y(test_idx, Eigen::all);

    MatrixType x_train = subset(x, train_idx);
    Eigen::MatrixXd y_train = y(train_idx, Eigen::all);

    return { x_train, y_train, x_test, y_test };
  }

  /**
   * @brief Get the number of folds
   *
   * @return size_t Number of folds
   */
  size_t numFolds() const { return n_folds; }

private:
  std::vector<std::vector<int>>
    folds;             ///< Vector of vectors containing indices for each fold
  std::size_t n_folds; ///< Number of folds

  /**
   * @brief Create random folds for cross-validation
   *
   * @param n_samples Total number of samples in the dataset
   * @param n_folds Number of folds to create
   * @param seed Random seed for reproducibility
   * @return std::vector<std::vector<int>> Vector of vectors containing indices
   * for each fold
   *
   * @internal
   * Implements stratified fold assignment by shuffling indices and
   * distributing them evenly across folds.
   * @endinternal
   */
  static std::vector<std::vector<int>> createFolds(int n_samples,
                                                   int n_folds,
                                                   uint64_t seed);
};

} // namespace slope
