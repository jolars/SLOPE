#include <algorithm>
#include <cassert>
#include <limits>
#include <slope/constants.h>
#include <slope/math.h>
#include <slope/sorted_l1_norm.h>
#include <slope/utils.h>

namespace slope {

double
SortedL1Norm::eval(const Eigen::VectorXd& beta,
                   const Eigen::ArrayXd& lambda) const
{
  assert(lambda.size() == beta.size() &&
         "Coefficient and lambda sizes must agree");

  Eigen::ArrayXd beta_abs = beta.array().abs();
  sort(beta_abs, true);
  return (beta_abs * lambda).sum();
}

Eigen::MatrixXd
SortedL1Norm::prox(const Eigen::VectorXd& beta,
                   const Eigen::ArrayXd& lambda) const
{
  // TODO: Avoid copying beta
  using namespace Eigen;

  assert(lambda.size() == beta.size() &&
         "Coefficient and lambda sizes must agree");

  ArrayXd beta_sign = beta.array().sign();
  VectorXd beta_copy = beta.array().abs();

  auto ord = sortIndex(beta_copy, true);
  permute(beta_copy, ord);

  int p = beta_copy.size();

  VectorXd s(p);
  VectorXd w(p);
  VectorXi idx_i(p);
  VectorXi idx_j(p);

  int k = 0;

  for (int i = 0; i < p; i++) {
    idx_i(k) = i;
    idx_j(k) = i;
    s(k) = beta_copy(i) - lambda(i);
    w(k) = s(k);

    while ((k > 0) && (w(k - 1) <= w(k))) {
      k--;
      idx_j(k) = i;
      s(k) += s(k + 1);
      w(k) = s(k) / (i - idx_i(k) + 1.0);
    }
    k++;
  }

  for (int j = 0; j < k; j++) {
    double d = std::max(w(j), 0.0);
    for (int i = idx_i(j); i <= idx_j(j); i++) {
      beta_copy(i) = d;
    }
  }

  // return order and sigsn
  inversePermute(beta_copy, ord);
  beta_copy.array() *= beta_sign;

  return beta_copy;
}

double
SortedL1Norm::dualNorm(const Eigen::VectorXd& gradient,
                       const Eigen::ArrayXd& lambda) const
{
  return dualNorm(gradient, lambda, 0.0);
}

double
SortedL1Norm::dualNorm(const Eigen::VectorXd& gradient,
                       const Eigen::ArrayXd& lambda,
                       const double denominator_floor) const
{
  Eigen::ArrayXd abs_gradient = gradient.cwiseAbs();
  sort(abs_gradient, true);

  assert(lambda.size() == gradient.size() &&
         "Gradient and lambda sizes must agree");

  const Eigen::ArrayXd cumulative_gradient = cumSum(abs_gradient);
  const Eigen::ArrayXd cumulative_lambda = cumSum(lambda);
  double result = 0.0;

  for (Eigen::Index i = 0; i < gradient.size(); ++i) {
    const double denominator =
      std::max(cumulative_lambda(i), denominator_floor);
    if (denominator == 0.0) {
      if (cumulative_gradient(i) > constants::MAX_DIV) {
        return std::numeric_limits<double>::infinity();
      }
      continue;
    }
    result = std::max(result, cumulative_gradient(i) / denominator);
  }

  return result;
}

} // namspace slope
