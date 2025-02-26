#include "regularization_sequence.h"
#include "qnorm.h"
#include "utils.h"
#include <cassert>

namespace slope {

Eigen::ArrayXd
lambdaSequence(const int p,
               const double q,
               const std::string& type,
               const int n,
               const double theta1,
               const double theta2)
{
  Eigen::ArrayXd lambda(p);

  validateOption(type, { "bh", "gaussian", "oscar", "lasso" }, "type");

  if (type == "bh" || type == "gaussian") {
    if (q <= 0 || q >= 1) {
      throw std::invalid_argument("q must be between 0 and 1");
    }

    for (int j = 0; j < p; ++j) {
      lambda(j) = normalQuantile(1.0 - (j + 1.0) * q / (2.0 * p));
    }

    if (type == "gaussian" && p > 1) {
      if (n <= 0) {
        throw std::invalid_argument(
          "n must be provided (and be positive) for type 'gaussian'");
      }

      double sum_sq = 0.0;

      for (int i = 1; i < p; ++i) {
        sum_sq += std::pow(lambda(i - 1), 2);
        double w = 1.0 / std::max(1.0, static_cast<double>(n - i - 1.0));

        lambda(i) *= std::sqrt(1.0 + w * sum_sq);
      }

      // Ensure non-increasing lambda
      for (int i = 1; i < p; ++i) {
        if (lambda(i) > lambda(i - 1)) {
          lambda(i) = lambda(i - 1);
        }
      }
    }
  } else if (type == "oscar") {
    if (theta1 < 0) {
      throw std::invalid_argument("theta1 must be non-negative");
    }
    if (theta2 < 0) {
      throw std::invalid_argument("theta2 must be non-negative");
    }
    lambda = theta1 + theta2 * (p - Eigen::ArrayXd::LinSpaced(p, 1, p));
  } else if (type == "lasso") {
    lambda.setOnes();
  }

  assert(lambda.minCoeff() > 0 && "lambda must be positive");
  assert(!lambda.isNaN().any() && "lambda must be finite");
  assert(lambda.size() == p && "lambda sequence is of right size");

  return lambda;
}

std::tuple<Eigen::ArrayXd, double, int>
regularizationPath(const Eigen::ArrayXd& alpha_in,
                   const Eigen::MatrixXd& gradient,
                   const SortedL1Norm& penalty,
                   const Eigen::ArrayXd& lambda,
                   const int n,
                   const int path_length,
                   double alpha_min_ratio)
{
  const int p = gradient.rows();
  const int m = gradient.cols();

  double alpha_max = penalty.dualNorm(gradient.reshaped(), lambda);

  if (alpha_in.size() != 0) {
    // User-supplied alpha sequence; just check it
    if (alpha_in.minCoeff() < 0) {
      throw std::invalid_argument("alpha must be non-negative");
    }
    if (!alpha_in.isFinite().all()) {
      throw std::invalid_argument("alpha must be finite");
    }

    return { alpha_in, alpha_max, alpha_in.size() };
  }

  if (alpha_min_ratio < 0) {
    alpha_min_ratio = n > p * m ? 1e-4 : 1e-2;
  }

  Eigen::ArrayXd alpha(path_length);

  double div = path_length - 1;

  for (int i = 0; i < path_length; ++i) {
    alpha(i) = alpha_max * std::pow(alpha_min_ratio, i / div);
  }

  return { alpha, alpha_max, path_length };
}

} // namespace slope
