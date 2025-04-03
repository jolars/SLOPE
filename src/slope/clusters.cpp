#include "clusters.h"
#include "math.h"
#include "utils.h"

namespace slope {

Clusters::Clusters(const Eigen::VectorXd& beta)
{
  update(beta);
}

bool
Clusters::hasAllZeros() const
{
  return !c.empty() && c[0] == 0 && c.size() == 1 && c_ind.empty();
}

std::vector<int>&
Clusters::getZeroIndices() const
{
  if (!zero_indices_valid) {
    zero_indices.clear();

    // Find indices not in c_ind (set difference)
    for (int j = 0; j < p; ++j) {
      if (std::find(c_ind.begin(), c_ind.end(), j) == c_ind.end()) {
        zero_indices.push_back(j);
      }
    }

    zero_indices_valid = true;
  }

  return zero_indices;
}

std::vector<int>::const_iterator
Clusters::cbegin(const int i) const
{
  if (i < static_cast<int>(c.size())) {
    return c_ind.cbegin() + c_ptr[i];
  }

  return getZeroIndices().cbegin();
}

std::vector<int>::const_iterator
Clusters::cend(const int i) const
{
  if (i < static_cast<int>(c.size())) {
    return c_ind.cbegin() + c_ptr[i + 1];
  }

  return getZeroIndices().cend();
}

std::vector<int>::iterator
Clusters::begin(const int i)
{
  if (i < static_cast<int>(c.size())) {
    return c_ind.begin() + c_ptr[i];
  }

  return getZeroIndices().begin();
}

std::vector<int>::iterator
Clusters::end(const int i)
{
  if (i < static_cast<int>(c.size())) {
    return c_ind.begin() + c_ptr[i + 1];
  }

  return getZeroIndices().end();
}

int
Clusters::cluster_size(const int i) const
{
  if (i < static_cast<int>(c.size())) {
    return c_ptr[i + 1] - c_ptr[i];
  }

  // For zero cluster, compute the size based on indices not in other clusters
  if (i == static_cast<int>(c.size()) && p > static_cast<int>(c_ind.size())) {
    return p - static_cast<int>(c_ind.size());
  }

  return 0; // Invalid cluster index
}

int
Clusters::pointer(const int i) const
{
  if (i < static_cast<int>(c_ptr.size())) {
    return c_ptr[i];
  }

  return c_ind.size(); // Return end of indices for virtual zero cluster
}

int
Clusters::n_clusters() const
{
  // Special case: if we have a coefficient of 0 in c, that means all values are
  // 0
  if (hasAllZeros()) {
    return 1;
  }

  // Regular case: count non-zero clusters, and add 1 if there are zero
  // coefficients
  return c.size() + (p > static_cast<int>(c_ind.size()) ? 1 : 0);
}

double
Clusters::coeff(const int i) const
{
  if (i < static_cast<int>(c.size())) {
    return c[i];
  }

  return 0.0; // Return 0 for the virtual zero cluster
}

void
Clusters::setCoeff(const int i, const double x)
{
  c[i] = x;
}

std::vector<double>
Clusters::coeffs() const
{
  std::vector<double> result = c;

  // Special case: if we have a coefficient of 0 in c and all values are zeros
  // don't add an additional zero
  if (hasAllZeros()) {
    return result; // Just return {0.0}
  }

  // Add zero coefficient if there's a zero cluster
  if (p > static_cast<int>(c_ind.size())) {
    result.push_back(0.0);
  }

  return result;
}

std::vector<int>
Clusters::indices() const
{
  return c_ind;
}

std::vector<int>
Clusters::pointers() const
{
  return c_ptr;
}

void
Clusters::update(const int old_index, const int new_index, const double c_new)
{
  // Invalidate zero indices cache - the cluster structure is changing
  zero_indices_valid = false;

  auto c_old = coeff(old_index);

  if (c_new != c_old) {
    if (c_new == coeff(new_index)) {
      merge(old_index, new_index);
    } else {
      setCoeff(old_index, c_new);
      if (old_index != new_index) {
        reorder(old_index, new_index);
      }
    }
  }
}

void
Clusters::update(const Eigen::VectorXd& beta)
{
  using sort_pair = std::pair<double, int>;

  p = beta.size();

  c.clear();
  c_ind.clear();
  c_ptr.clear();
  signs.clear();

  // Invalidate zero indices cache
  zero_indices_valid = false;

  std::vector<sort_pair> sorted;
  sorted.reserve(p);
  signs.reserve(p);

  for (int i = 0; i < beta.size(); ++i) {
    signs.emplace_back(sign(beta(i)));
    double abs_val = std::abs(beta(i));
    if (abs_val > 0) {
      sorted.emplace_back(abs_val, i);
    }
  }

  // Special case: all values are zero
  if (sorted.empty() && p > 0) {
    c.push_back(0.0); // Add a single zero coefficient
    c_ptr.push_back(0);
    c_ptr.push_back(0);

    return;
  }

  std::sort(sorted.begin(), sorted.end(), std::greater<sort_pair>());

  c_ind.reserve(sorted.size());
  for (const auto& sorted_i : sorted) {
    c_ind.emplace_back(sorted_i.second);
  }

  // Extract unique coefficients while preserving order
  std::vector<sort_pair> sorted_unique;
  sorted_unique.reserve(sorted.size()); // At most, all values could be unique

  for (auto it = sorted.begin(); it != sorted.end();) {
    const double current_value = it->first;
    sorted_unique.emplace_back(*it);
    it = std::find_if(it, sorted.end(), [current_value](const sort_pair& elem) {
      return elem.first != current_value;
    });
  }

  c.reserve(sorted_unique.size());
  for (const auto& sorted_unique_i : sorted_unique) {
    c.emplace_back(sorted_unique_i.first);
  }

  c_ptr.reserve(c.size() + 1);
  c_ptr.emplace_back(0);

  auto range_start = sorted.begin();
  for (const auto& c_i : c) {
    auto range_end =
      std::find_if(range_start, sorted.end(), [&c_i](const sort_pair& x) {
        return x.first != c_i;
      });
    c_ptr.emplace_back(std::distance(sorted.begin(), range_end));
    range_start = range_end;
  }
}

void
Clusters::reorder(const int old_index, const int new_index)
{
  auto c_size = cluster_size(old_index);

  // update coefficients
  move_elements(c, old_index, new_index, 1);

  // update indices
  move_elements(c_ind, pointer(old_index), pointer(new_index), c_size);

  // update pointers
  if (new_index < old_index) {
    move_elements(c_ptr, old_index + 1, new_index + 1, 1);

    std::for_each(c_ptr.begin() + new_index + 1,
                  c_ptr.begin() + old_index + 2,
                  [c_size](int& x) { x += c_size; });

    c_ptr[new_index + 1] = c_ptr[new_index] + c_size;
  } else {
    move_elements(c_ptr, old_index, new_index, 1);

    std::for_each(c_ptr.begin() + old_index,
                  c_ptr.begin() + new_index,
                  [c_size](int& x) { x -= c_size; });
    c_ptr[new_index] = c_ptr[new_index + 1] - c_size;
  }
}

void
Clusters::merge(const int old_index, const int new_index)
{
  auto c_size = cluster_size(old_index);

  // update coefficients
  c.erase(c.cbegin() + old_index);

  // update indices
  move_elements(c_ind, pointer(old_index), pointer(new_index), c_size);

  // update pointers
  if (new_index < old_index) {
    std::for_each(c_ptr.begin() + new_index + 1,
                  c_ptr.begin() + old_index + 1,
                  [c_size](int& x) { x += c_size; });
  } else {
    std::for_each(c_ptr.begin() + old_index + 1,
                  c_ptr.begin() + new_index + 1,
                  [c_size](int& x) { x -= c_size; });
  }

  c_ptr.erase(c_ptr.begin() + old_index + 1);
}

std::vector<std::vector<int>>
Clusters::getClusters() const
{
  std::vector<std::vector<int>> clusters;
  clusters.reserve(n_clusters());

  // Special case for all zeros vector
  if (hasAllZeros()) {
    std::vector<int> zero_cluster;
    for (int i = 0; i < p; ++i) {
      zero_cluster.push_back(i);
    }
    clusters.push_back(zero_cluster);
    return clusters;
  }

  for (int i = 0; i < n_clusters(); ++i) {
    clusters.emplace_back(cbegin(i), cend(i));
  }

  return clusters;
}

Eigen::SparseMatrix<int>
Clusters::patternMatrix() const
{
  std::vector<Eigen::Triplet<int>> triplets;
  triplets.reserve(p);

  int n_cols = 0;

  for (int k = 0; k < n_clusters(); ++k) {
    if (coeff(k) == 0) {
      // Skip the zero cluster
      break;
    }

    n_cols++;

    for (auto it = cbegin(k); it != cend(k); ++it) {
      int ind = *it;
      int s = signs[ind];
      triplets.emplace_back(ind, k, s);
    }
  }

  Eigen::SparseMatrix<int> out(p, n_cols);
  out.setFromTriplets(triplets.begin(), triplets.end());

  return out;
}

Eigen::SparseMatrix<int>
patternMatrix(const Eigen::VectorXd& beta)
{
  Clusters clusters(beta);
  return clusters.patternMatrix();
}

} // namespace slope
