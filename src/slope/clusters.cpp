#include "clusters.h"
#include "utils.h"

namespace slope {

Clusters::Clusters(const Eigen::VectorXd& beta)
{
  update(beta);
}

std::vector<int>::iterator
Clusters::begin(const int i)
{
  return c_ind.begin() + c_ptr[i];
}

std::vector<int>::iterator
Clusters::end(const int i)
{
  return c_ind.begin() + c_ptr[i + 1];
}

std::vector<int>::const_iterator
Clusters::cbegin(const int i) const
{
  return c_ind.cbegin() + c_ptr[i];
}

std::vector<int>::const_iterator
Clusters::cend(const int i) const
{
  return c_ind.cbegin() + c_ptr[i + 1];
}

int
Clusters::cluster_size(const int i) const
{
  return c_ptr[i + 1] - c_ptr[i];
}

int
Clusters::pointer(const int i) const
{
  return c_ptr[i];
}

int
Clusters::n_clusters() const
{
  return c.size();
}

double
Clusters::coeff(const int i) const
{
  return c[i];
}

void
Clusters::setCoeff(const int i, const double x)
{
  c[i] = x;
}

std::vector<double>
Clusters::coeffs() const
{
  return c;
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

  c.clear();
  c_ind.clear();
  c_ptr.clear();

  std::vector<sort_pair> sorted;
  sorted.reserve(beta.size());

  for (int i = 0; i < beta.size(); ++i) {
    sorted.emplace_back(std::abs(beta(i)), i);
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

  for (int i = 0; i < n_clusters(); ++i) {
    clusters.emplace_back(cbegin(i), cend(i));
  }

  return clusters;
}

} // namespace slope
