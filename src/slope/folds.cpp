#include "folds.h"
#include <algorithm>
#include <random>

namespace slope {

const std::vector<int>&
Folds::getTestIndices(size_t fold_idx, size_t rep_idx) const
{
  return folds[rep_idx][fold_idx];
}

std::vector<int>
Folds::getTrainingIndices(size_t fold_idx, size_t rep_idx) const
{
  std::vector<int> train_indices;
  for (size_t i = 0; i < n_folds; ++i) {
    if (i != fold_idx) {
      const auto& fold = folds[rep_idx][i];
      train_indices.insert(train_indices.end(), fold.begin(), fold.end());
    }
  }
  return train_indices;
}

std::vector<std::vector<int>>
Folds::createFolds(int n, int n_folds, uint64_t random_seed)
{
  // Initialize random number generator
  std::mt19937 generator(random_seed);

  // Create and shuffle indices
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), generator);

  // Create folds
  std::vector<std::vector<int>> folds(n_folds);

  // Calculate base fold size and remainder
  int base_fold_size = n / n_folds;
  int remainder = n % n_folds;

  // Current position in indices
  int current_pos = 0;

  // Distribute indices across folds
  for (int fold = 0; fold < n_folds; ++fold) {
    // Add one extra element to early folds if we have remainder
    int fold_size = base_fold_size + (fold < remainder ? 1 : 0);

    // Fill this fold
    folds[fold].reserve(fold_size);
    for (int i = 0; i < fold_size; ++i) {
      folds[fold].push_back(indices[current_pos++]);
    }
  }

  return folds;
}

} // namespace slope
