#include "screening.h"
#include "slope/utils.h"
#include <Eigen/Core>

namespace slope {

typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;
typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXXb;

std::vector<int>
activeSet(const Eigen::MatrixXd& beta)
{
  ArrayXb active = (beta.array() != 0.0).rowwise().any();

  return which(active);
}

std::vector<int>
strongSet(const Eigen::MatrixXd& gradient_prev,
          const Eigen::ArrayXd& lambda,
          const Eigen::ArrayXd& lambda_prev)
{
  using Eigen::VectorXd;
  using Eigen::VectorXi;

  int p = gradient_prev.rows();
  int m = gradient_prev.cols();

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

  while (i + k < p * m) {
    s += tmp(k + i);

    if (s >= 0) {
      k = k + i + 1;
      i = 0;
      s = 0;
    } else {
      i++;
    }
  }

  ArrayXXb active_set = ArrayXXb::Zero(p * m, 1);

  active_set.topRows(k).setOnes();

  // reset order
  active_set(ord, 0) = active_set.col(0).eval();

  active_set.resize(p, m);

  // TODO: Don't use rowwise here. Coefficients should be screened elementwise.
  ArrayXb active = active_set.array().rowwise().any();

  return which(active);
}

} // namespace slope
