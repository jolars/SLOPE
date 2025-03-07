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

  assert(out == out && "Loss is NaN");

  out -= (y.array() * eta.array()).sum() / n;

  return out;
}

double
Multinomial::dual(const Eigen::MatrixXd& theta,
                  const Eigen::MatrixXd& y,
                  const Eigen::VectorXd&)
{
  const Eigen::MatrixXd r = theta + y;

  return -(r.array() * r.array().max(constants::P_MIN).log()).sum() / y.rows();
}

Eigen::MatrixXd
Multinomial::residual(const Eigen::MatrixXd& eta, const Eigen::MatrixXd& y)
{
  return softmax(eta) - y;
}

Eigen::MatrixXd
Multinomial::preprocessResponse(const Eigen::MatrixXd& y)
{
  const int n = y.rows();
  int m = y.cols();

  Eigen::MatrixXd result;

  if (m == 1) {
    // Y is a column vector, expect integers representing classes
    m = y.array().maxCoeff() + 1; // Assuming classes are 0-based

    if (m == 1) {
      throw std::invalid_argument("Only one class found in response");
    }

    result = Eigen::MatrixXd::Zero(n, m);

    for (int i = 0; i < n; i++) {
      int class_label = static_cast<int>(y(i, 0));
      if (class_label < 0) {
        throw std::invalid_argument(
          "Class labels must be consecutive integers starting from 0");
      }

      result(i, class_label) = 1.0;
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

    result = y;
  }

  return result;
}

void
Multinomial::updateWeightsAndWorkingResponse(Eigen::VectorXd&,
                                             Eigen::VectorXd&,
                                             const Eigen::VectorXd&,
                                             const Eigen::VectorXd&)
{
  throw std::runtime_error("Multinomial loss does not currently support IRLS");
}

Eigen::MatrixXd
Multinomial::link(const Eigen::MatrixXd& mu)
{
  return mu.unaryExpr([](const double& x) {
    return logit(std::clamp(x, constants::P_MIN, constants::P_MAX));
  });
}

Eigen::MatrixXd
Multinomial::inverseLink(const Eigen::MatrixXd& eta)
{
  return softmax(eta);
}

Eigen::MatrixXd
Multinomial::predict(const Eigen::MatrixXd& eta)
{
  Eigen::MatrixXd prob = inverseLink(eta);

  // Find the class with the highest probability
  Eigen::VectorXd out(eta.rows());

  for (int i = 0; i < eta.rows(); i++) {
    out(i) = whichMax(prob.row(i));
  }

  return out;
}

// TODO: Consider adjusting the coefficients somehow.

} // namespace slope
