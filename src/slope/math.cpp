#include "math.h"
#include "constants.h"

namespace slope {

Eigen::VectorXd
logSumExp(const Eigen::MatrixXd& a)
{
  Eigen::ArrayXd max_vals = a.rowwise().maxCoeff();
  Eigen::ArrayXd sum_exp = (-max_vals).min(constants::MAX_EXP).exp() +
                           (a.colwise() - max_vals.matrix())
                             .array()
                             .min(constants::MAX_EXP)
                             .exp()
                             .rowwise()
                             .sum();

  return max_vals + sum_exp.max(constants::P_MIN).log();
}

Eigen::MatrixXd
softmax(const Eigen::MatrixXd& eta)
{
  return (eta.colwise() - logSumExp(eta))
    .array()
    .min(constants::MAX_EXP)
    .exp()
    .max(constants::P_MIN)
    .min(constants::P_MAX);
}

std::vector<int>
setUnion(const std::vector<int>& a, const std::vector<int>& b)
{
  assert(std::is_sorted(a.begin(), a.end()) &&
         "First argument to setUnion must be sorted");
  assert(std::is_sorted(b.begin(), b.end()) &&
         "Second argument to setUnion must be sorted");

  std::vector<int> out;
  std::set_union(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  return out;
}

std::vector<int>
setDiff(const std::vector<int>& a, const std::vector<int>& b)
{
  assert(std::is_sorted(a.begin(), a.end()) &&
         "First argument to setDiff must be sorted");
  assert(std::is_sorted(b.begin(), b.end()) &&
         "Second argument to setDiff must be sorted");

  std::vector<int> out;
  std::set_difference(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  return out;
}

Eigen::ArrayXd
geomSpace(const double start, const double end, const int n)
{
  if (n == 1) {
    return Eigen::ArrayXd::Constant(1, start);
  } else {
    return Eigen::ArrayXd::LinSpaced(n, std::log(start), std::log(end)).exp();
  }
}

} // namespace slope
