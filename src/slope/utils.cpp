#include "utils.h"
#include <stdexcept>

namespace slope {

void
validateOption(const std::string& value,
               const std::set<std::string>& valid_options,
               const std::string& parameter_name)
{
  if (valid_options.find(value) == valid_options.end()) {
    std::string valid_list =
      std::accumulate(std::next(valid_options.begin()),
                      valid_options.end(),
                      std::string("'") + *valid_options.begin() + "'",
                      [](const std::string& a, const std::string& b) {
                        return a + ", '" + b + "'";
                      });

    throw std::invalid_argument("Invalid " + parameter_name + ": '" + value +
                                "'. Must be one of: " + valid_list);
  }
}

Eigen::MatrixXd
subset(const Eigen::MatrixXd& x, const std::vector<int>& indices)
{
  return x(indices, Eigen::all);
}

Eigen::SparseMatrix<double>
subset(const Eigen::SparseMatrix<double>& x, const std::vector<int>& indices)
{
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(x.nonZeros());

  for (int j = 0; j < x.cols(); ++j) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(x, j); it; ++it) {
      auto it_idx = std::find(indices.begin(), indices.end(), it.row());

      if (it_idx != indices.end()) {
        int new_row = std::distance(indices.begin(), it_idx);
        triplets.emplace_back(new_row, j, it.value());
      }
    }
  }
  Eigen::SparseMatrix<double> out(indices.size(), x.cols());
  out.setFromTriplets(triplets.begin(), triplets.end());

  return out;
}

} // namespace slope
