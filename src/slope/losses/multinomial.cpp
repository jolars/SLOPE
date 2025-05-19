#include "multinomial.h"
#include "../constants.h"
#include "../math.h"
#include "../utils.h"
#include <unordered_set>

namespace slope {

double
Multinomial::loss(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  int n = y.rows();

  double out = logSumExp(eta).mean();

  assert(eta.allFinite());
  assert(out == out && "Loss is NaN");

  out -= (y.array() * eta.array()).sum() / n;

  return out;
}

double
Multinomial::dual(const Eigen::MatrixXd& theta,
                  const Eigen::MatrixXd& y,
                  const Eigen::VectorXd&)
{
  int n = y.rows();
  int p = y.cols();

  Eigen::ArrayXXd eta = link(theta + y);

  // TODO: Find out if this formulation can be improved numerically
  double out =
    logSumExp(eta).mean() - (eta * (y.array() + theta.array())).sum() / n;

  return out;
}

Eigen::MatrixXd
Multinomial::preprocessResponse(const Eigen::MatrixXd& y)
{
  const int n = y.rows();

  Eigen::MatrixXd result;

  if (y.cols() == 1) {
    // y is a column vector, expect integers representing classes
    int m = y.array().maxCoeff(); // Assuming classes are 0-based

    if (m == 0) {
      throw std::invalid_argument("Only one class found in response");
    }

    result = Eigen::MatrixXd::Zero(n, m);

    for (int i = 0; i < n; i++) {
      int class_label = static_cast<int>(y(i, 0));
      if (class_label < 0) {
        throw std::invalid_argument(
          "Class labels must be consecutive integers starting from 0");
      }

      if (class_label < m) {
        result(i, class_label) = 1.0;
      }
    }
  } else {
    // Y is a matrix, expect one-hot encoding
    auto y_unique = unique(y);

    if (y_unique.size() > 2) {
      throw std::invalid_argument(
        "Expected binary labels (0/1) but found more than two unique values");
    }

    for (const auto& val : y_unique) {
      if (val != 0.0 && val != 1.0) {
        throw std::invalid_argument(
          "Expected binary labels with values 0 and 1 only");
      }
    }

    return y;
  }

  return result;
}

Eigen::MatrixXd
Multinomial::hessianDiagonal(const Eigen::MatrixXd& eta)
{
  Eigen::MatrixXd pr = inverseLink(eta);

  return pr.array() * (1.0 - pr.array());
}

Eigen::MatrixXd
Multinomial::link(const Eigen::MatrixXd& mu)
{
  Eigen::MatrixXd out(mu.rows(), mu.cols());

  for (int i = 0; i < mu.rows(); i++) {
    double row_sum = mu.row(i).sum();
    double ref_val = 1.0 - row_sum;

    for (int j = 0; j < mu.cols(); j++) {
      double numerator =
        std::clamp(mu(i, j), constants::P_MIN, constants::P_MAX);
      double denominator =
        std::clamp(ref_val, constants::P_MIN, constants::P_MAX);
      out(i, j) = std::log(numerator) - std::log(denominator);
    }
  }

  return out;
}

Eigen::MatrixXd
Multinomial::inverseLink(const Eigen::MatrixXd& eta)
{
  return softmax(eta);
}

Eigen::MatrixXd
Multinomial::predict(const Eigen::MatrixXd& eta)
{
  int n = eta.rows();
  int m = eta.cols();

  // Directly compute probabilities with the last column as reference class
  Eigen::MatrixXd prob = softmax(eta);

  // Find the class with the highest probability
  Eigen::VectorXd out(n);
  for (int i = 0; i < n; i++) {
    double sum_prob = prob.row(i).sum();
    double best_prob = prob.row(i).maxCoeff();
    double ref_prob = 1.0 - sum_prob;

    out(i) = best_prob > ref_prob ? whichMax(prob.row(i)) : m;
  }

  return out;
}

// TODO: Consider adjusting the coefficients somehow.

} // namespace slope
