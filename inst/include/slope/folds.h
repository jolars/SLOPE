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
   * @brief Default constructor
   *
   * Creates an empty Folds object that can be assigned to later.
   */
  Folds()
    : n_folds(0)
    , n_repeats(0)
  {
  }

  /**
   * @brief Constructor for generating random folds with optional repetitions
   *
   * @param n_samples Total number of samples in the dataset
   * @param n_folds Number of folds to create
   * @param n_repeats Number of repetitions (default: 1)
   * @param seed Random seed for reproducibility
   *
   * Creates k-fold partitions by randomly assigning samples to folds.
   * When n_repeats > 1, multiple sets of folds are created with different
   * random seeds.
   */
  Folds(int n_samples, int n_folds, int n_repeats = 1, uint64_t seed = 42)
    : n_folds(n_folds)
    , n_repeats(n_repeats)
  {
    folds.resize(n_repeats);
    for (int rep = 0; rep < n_repeats; ++rep) {
      // Use a different seed for each repetition
      uint64_t rep_seed = seed + rep;
      folds[rep] = createFolds(n_samples, n_folds, rep_seed);
    }
  }

  /**
   * @brief Constructor for user-provided folds
   *
   * @param user_folds Vector of vectors containing indices for each fold
   *
   * Uses the provided fold assignments instead of creating random folds.
   */
  explicit Folds(const std::vector<std::vector<int>>& user_folds)
    : n_folds(user_folds.size())
    , n_repeats(1)
  {
    folds.emplace_back(user_folds);
  }

  /**
   * @brief Constructor for user-provided repeated folds
   *
   * @param user_folds Vector of fold sets for multiple repetitions
   */
  explicit Folds(const std::vector<std::vector<std::vector<int>>>& user_folds)
    : folds(user_folds)
    , n_folds(user_folds[0].size())
    , n_repeats(user_folds.size())
  {
  }

  /**
   * @brief Get test indices for a specific fold and repetition
   *
   * @param fold_idx Index of the fold to use as test set
   * @param rep_idx Index of the repetition (default: 0)
   * @return const std::vector<int>& Reference to the vector of test indices
   */
  const std::vector<int>& getTestIndices(size_t fold_idx,
                                         size_t rep_idx = 0) const;

  /**
   * @brief Get training indices for a specific fold and repetition
   *
   * @param fold_idx Index of the fold to exclude from training
   * @param rep_idx Index of the repetition (default: 0)
   * @return std::vector<int> Vector containing all indices except those in the
   * specified fold
   */
  std::vector<int> getTrainingIndices(size_t fold_idx,
                                      size_t rep_idx = 0) const;

  /**
   * @brief Split data into training and test sets for a specific fold and
   * repetition
   *
   * @tparam T Type of design matrix (supports both dense and sparse matrices)
   * @param x Input feature matrix
   * @param y Response matrix
   * @param fold_idx Index of the fold to use as test set
   * @param rep_idx Index of the repetition (default: 0)
   * @return std::tuple containing (x_train, y_train, x_test, y_test)
   *
   * This method creates training and test datasets by subsetting the original
   * data according to the specified fold indices. For dense matrices, it
   * creates copies of the data. For sparse matrices, it creates efficiently
   * constructed sparse matrix subsets.
   */
  template<typename T>
  auto split(Eigen::EigenBase<T>& x,
             const Eigen::MatrixXd& y,
             size_t fold_idx,
             size_t rep_idx = 0) const
  {
    auto test_idx = getTestIndices(fold_idx, rep_idx);
    auto train_idx = getTrainingIndices(fold_idx, rep_idx);

    auto x_test = subset(x.derived(), test_idx);
    Eigen::MatrixXd y_test = y(test_idx, Eigen::all);

    auto x_train = subset(x.derived(), train_idx);
    Eigen::MatrixXd y_train = y(train_idx, Eigen::all);

    return std::make_tuple(x_train, y_train, x_test, y_test);
  }

  /**
   * @brief Get the number of folds
   *
   * @return size_t Number of folds
   */
  size_t numFolds() const { return n_folds; }

  /**
   * @brief Get the number of repetitions
   *
   * @return size_t Number of repetitions
   */
  size_t numRepetitions() const { return n_repeats; }

  /**
   * @brief Get the total number of folds (repetitions * folds)
   *
   * @return size_t Number of evaluations
   */
  size_t numEvals() const { return n_repeats * n_folds; }

  std::vector<std::vector<std::vector<int>>>
    folds;               ///< Indices for each fold in each repetition
  std::size_t n_folds;   ///< Number of folds
  std::size_t n_repeats; ///< Number of repetitions

private:
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
