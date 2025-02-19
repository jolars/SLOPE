#include "math.h"
#include "constants.h"

namespace slope {

Eigen::VectorXd
logSumExp(const Eigen::MatrixXd& a)
{
  Eigen::VectorXd max_vals = a.rowwise().maxCoeff();
  Eigen::ArrayXd sum_exp =
    (a.colwise() - max_vals).array().exp().rowwise().sum();

  return max_vals.array() + sum_exp.max(constants::P_MIN).log();
}

Eigen::MatrixXd
softmax(const Eigen::MatrixXd& a)
{
  Eigen::VectorXd shift = a.rowwise().maxCoeff();
  Eigen::MatrixXd exp_a = (a.colwise() - shift).array().exp();
  Eigen::ArrayXd row_sums = exp_a.rowwise().sum();
  return exp_a.array().colwise() / row_sums;
}

std::vector<int>
setUnion(const std::vector<int>& a, const std::vector<int>& b)
{
  std::vector<int> out;
  std::set_union(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  return out;
}

std::vector<int>
setDiff(const std::vector<int>& a, const std::vector<int>& b)
{
  std::vector<int> out;
  std::set_difference(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  return out;
}

} // namespace slope
