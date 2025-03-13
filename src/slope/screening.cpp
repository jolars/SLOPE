#include "screening.h"
#include "utils.h"
#include <Eigen/Core>

namespace slope {

typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;
typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXXb;

std::vector<int>
activeSet(const Eigen::VectorXd& beta)
{
  ArrayXb active = beta.array() != 0.0;

  return which(active);
}

std::vector<int>
strongSet(const Eigen::VectorXd& gradient_prev,
          const Eigen::ArrayXd& lambda,
          const Eigen::ArrayXd& lambda_prev)
{
  using Eigen::VectorXd;
  using Eigen::VectorXi;

  int pm = gradient_prev.size();

  assert(lambda_prev.size() == lambda.size() &&
         "lambda_prev and lambda must have the same length");
  assert((lambda <= lambda_prev).all() &&
         "New lambda values must be smaller than or equal to previous values");

  const VectorXd abs_grad = gradient_prev.reshaped().cwiseAbs();
  std::vector<int> ord = sortIndex(abs_grad, true);

  assert(abs_grad.size() == lambda.size());

  const VectorXd tmp =
    abs_grad(ord).array().eval() + lambda_prev - 2.0 * lambda;

  int i = 0;
  int k = 0;

  double s = 0;

  while (i + k < pm) {
    s += tmp(k + i);

    if (s >= 0) {
      k = k + i + 1;
      i = 0;
      s = 0;
    } else {
      i++;
    }
  }

  ArrayXb active_set = ArrayXb::Zero(pm);
  active_set.head(k).setOnes();

  // restore order
  inversePermute(active_set, ord);

  return which(active_set);
}

} // namespace slope
