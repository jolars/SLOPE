#include "clusters.h"
#include "math.h"
#include "utils.h"

namespace slope {

Clusters::Clusters(const Eigen::VectorXd& beta)
{
  update(beta);
}

std::vector<int>::const_iterator
Clusters::cbegin(const int i) const
{
  assert(i >= 0 && i < size());
  return c_ind.cbegin() + this->pointer(i);
}

std::vector<int>::const_iterator
Clusters::cend(const int i) const
{
  assert(i >= 0 && i < size());
  return c_ind.cbegin() + this->pointer(i + 1);
}

std::vector<int>::iterator
Clusters::begin(const int i)
{
  assert(i >= 0 && i < size());
  return c_ind.begin() + this->pointer(i);
}

std::vector<int>::iterator
Clusters::end(const int i)
{
  assert(i >= 0 && i < size());
  return c_ind.begin() + this->pointer(i + 1);
}

int
Clusters::cluster_size(const int i) const
{
  assert(i >= 0 && i < size());
  return this->pointer(i + 1) - this->pointer(i);
}

int
Clusters::pointer(const int i) const
{
  assert(i >= 0 && i <= size());
  return c_ptr[i];
}

int
Clusters::size() const
{
  return c.size();
}

double
Clusters::coeff(const int i) const
{
  assert(i >= 0 && i < size());
  return c[i];
}

void
Clusters::setCoeff(const int i, const double x)
{
  assert(i >= 0 && i < size());
  c[i] = x;
}

const std::vector<double>&
Clusters::coeffs() const
{
  return c;
}

const std::vector<int>&
Clusters::indices() const
{
  return c_ind;
}

const std::vector<int>&
Clusters::pointers() const
{
  return c_ptr;
}

void
Clusters::update(const int old_index, const int new_index, const double c_new)
{
  assert(old_index < size());
  assert(new_index <= size());

  auto c_old = coeff(old_index);

  if (c_new != c_old) {
    if (c_new == 0) {
      // Get the size and pointer of the cluster to be removed
      auto cluster_size_val = cluster_size(old_index);

      // Remove indices in the cluster
      c_ind.erase(c_ind.begin() + pointer(old_index),
                  c_ind.begin() + pointer(old_index + 1));

      // Update pointers after the removed cluster
      c_ptr.erase(c_ptr.begin() + old_index + 1);

      for (size_t i = old_index + 1; i < c_ptr.size(); ++i) {
        c_ptr[i] -= cluster_size_val;
      }

      c.erase(c.begin() + old_index);
    } else if (c_new == coeff(new_index)) {
      merge(old_index, new_index);
    } else {
      setCoeff(old_index, c_new);
      if (old_index != new_index) {
        reorder(old_index, new_index);
      }
    }
  }

  assert(c_ptr.size() == c.size() + 1);
  assert(std::is_sorted(c.begin(), c.end(), std::greater<>{}));
}

void
Clusters::update(const Eigen::VectorXd& beta)
{
  using sort_pair = std::pair<double, int>;

  p = beta.size();

  c.clear();
  c_ind.clear();
  c_ptr.clear();

  std::vector<sort_pair> sorted;
  sorted.reserve(p);

  for (int i = 0; i < beta.size(); ++i) {
    double abs_val = std::abs(beta(i));
    if (abs_val > 0) {
      sorted.emplace_back(abs_val, i);
    }
  }

  sorted.shrink_to_fit();

  // Special case: all values are zero
  if (sorted.empty()) {
    c_ptr.push_back(0); // Initialize with a single pointer at 0
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

  assert(c_ptr.size() == c.size() + 1);
  assert(std::is_sorted(c.begin(), c.end(), std::greater<>{}));
}

void
Clusters::reorder(const int old_index, const int new_index)
{
  auto c_size = cluster_size(old_index);
  int old_ptr = pointer(old_index); // Save pointer for cluster being moved
  int new_ptr = pointer(new_index); // Save destination pointer

  // update coefficients
  move_elements(c, old_index, new_index, 1);

  assert(std::is_sorted(c.begin(), c.end(), std::greater<>{}));

  // update pointers
  if (new_index < old_index) {
    std::rotate(c_ind.begin() + pointer(new_index),
                c_ind.begin() + pointer(old_index),
                c_ind.begin() + pointer(old_index + 1));

    move_elements(c_ptr, old_index + 1, new_index + 1, 1);

    std::for_each(c_ptr.begin() + new_index + 1,
                  c_ptr.begin() + old_index + 2,
                  [c_size](int& x) { x += c_size; });

    c_ptr[new_index + 1] = c_ptr[new_index] + c_size;
  } else {
    std::rotate(c_ind.begin() + pointer(old_index),
                c_ind.begin() + pointer(old_index + 1),
                c_ind.begin() + pointer(new_index + 1));

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
  assert(old_index >= 0 && "old_index must be non-negative");
  assert(new_index >= 0 && "new_index must be non-negative");
  assert(old_index != new_index && "Cannot merge a cluster with itself");
  assert(old_index < size() && "old_index out of bounds");
  assert(new_index < size() && "new_index out of bounds");

  auto c_size_old = cluster_size(old_index);
  auto c_size_new = cluster_size(new_index);

  auto c_ind_begin = c_ind.begin();

  assert(pointer(old_index) >= 0 &&
         pointer(old_index) <= static_cast<int>(c_ind.size()));
  assert(pointer(old_index + 1) >= 0 &&
         pointer(old_index + 1) <= static_cast<int>(c_ind.size()));
  assert(pointer(new_index) >= 0 &&
         pointer(new_index) <= static_cast<int>(c_ind.size()));
  assert(pointer(new_index + 1) >= 0 &&
         pointer(new_index + 1) <= static_cast<int>(c_ind.size()));

  assert(c_ind_begin + pointer(old_index) >= c_ind.begin() &&
         c_ind_begin + pointer(old_index) <= c_ind.end());
  assert(c_ind_begin + pointer(old_index + 1) >= c_ind.begin() &&
         c_ind_begin + pointer(old_index + 1) <= c_ind.end());
  assert(c_ind_begin + pointer(new_index + 1) >= c_ind.begin() &&
         c_ind_begin + pointer(new_index + 1) <= c_ind.end());

  // update pointers
  if (new_index < old_index) {
    assert(pointer(new_index + 1) <= static_cast<int>(c_ind.size()));
    assert(pointer(old_index) < pointer(old_index + 1) &&
           "First two iterators must form a valid range");
    assert(pointer(old_index) <= pointer(old_index + 1) &&
           "Second iterator must be before or equal to third");
    assert(pointer(new_index + 1) <= static_cast<int>(c_ind.size()));

    std::rotate(c_ind_begin + pointer(new_index),
                c_ind_begin + pointer(old_index),
                c_ind_begin + pointer(old_index + 1));
    std::for_each(c_ptr.begin() + new_index + 1,
                  c_ptr.begin() + old_index + 1,
                  [c_size_old](int& x) { x += c_size_old; });
  } else {
    assert(pointer(old_index + 1) <= static_cast<int>(c_ind.size()));
    assert(pointer(old_index) < pointer(old_index + 1) &&
           "First two iterators must form a valid range");
    assert(pointer(old_index + 1) <= pointer(new_index + 1) &&
           "Second iterator must be before or equal to third");
    assert(pointer(old_index + 1) <= static_cast<int>(c_ind.size()));

    std::rotate(c_ind_begin + pointer(old_index),
                c_ind_begin + pointer(old_index + 1),
                c_ind_begin + pointer(new_index + 1));
    std::for_each(c_ptr.begin() + old_index + 1,
                  c_ptr.begin() + new_index + 1,
                  [c_size_old](int& x) { x -= c_size_old; });
  }

  c_ptr.erase(c_ptr.begin() + old_index + 1);

  // update coefficients
  c.erase(c.cbegin() + old_index);

  assert(c_ptr.size() == c.size() + 1 &&
         "Pointer array size mismatch after merge");
}

std::vector<std::vector<int>>
Clusters::getClusters() const
{
  std::vector<std::vector<int>> clusters;
  clusters.reserve(size());

  for (int i = 0; i < size(); ++i) {
    clusters.emplace_back(cbegin(i), cend(i));
  }

  return clusters;
}

// TODO: Template or overload for sparse matrices
Eigen::SparseMatrix<int>
patternMatrix(const Eigen::VectorXd& beta)
{
  Clusters clusters(beta);

  int p = beta.size();
  Eigen::SparseMatrix<int> out(p, clusters.size());

  if (clusters.size() == 0) {
    return out; // Return empty matrix if no clusters
  }

  std::vector<Eigen::Triplet<int>> triplets;
  triplets.reserve(p);

  for (int k = 0; k < clusters.size(); ++k) {
    for (auto it = clusters.cbegin(k); it != clusters.cend(k); ++it) {
      int ind = *it;
      int s = std::copysign(1.0, beta(ind));
      triplets.emplace_back(ind, k, s);
    }
  }

  out.setFromTriplets(triplets.begin(), triplets.end());

  return out;
}

} // namespace slope
