#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <slope/constants.h>
#include <slope/losses/multinomial.h>
#include <slope/math.h>
#include <slope/utils.h>
#include <stdexcept>
#include <unordered_set>

namespace {

double
xLogX(const double x)
{
  return x == 0.0 ? 0.0 : x * std::log(x);
}

} // namespace

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
  if (theta.rows() != y.rows() || theta.cols() != y.cols()) {
    throw std::invalid_argument(
      "Multinomial dual point and response dimensions must match");
  }

  const Eigen::ArrayXXd mean = theta.array() + y.array();
  double entropy = 0.0;

  for (Eigen::Index i = 0; i < mean.rows(); ++i) {
    double reference = 1.0;
    for (Eigen::Index k = 0; k < mean.cols(); ++k) {
      const double value = mean(i, k);
      if (!std::isfinite(value) || value < 0.0) {
        throw std::domain_error(
          "Multinomial dual rows must be finite and lie in the simplex");
      }
      entropy += xLogX(value);
      reference -= value;
    }
    if (!std::isfinite(reference) || reference < 0.0) {
      throw std::domain_error(
        "Multinomial dual rows must be finite and lie in the simplex");
    }
    entropy += xLogX(reference);
  }

  return -entropy / y.rows();
}

Eigen::MatrixXd
Multinomial::dualPoint(const Eigen::MatrixXd& eta,
                       const Eigen::MatrixXd& y,
                       const bool fit_intercept)
{
  Eigen::MatrixXd theta = residual(eta, y);
  if (!fit_intercept) {
    return theta;
  }
  if (eta.rows() != y.rows() || eta.cols() != y.cols()) {
    throw std::invalid_argument(
      "Multinomial dual point and response dimensions must match");
  }

  theta.rowwise() -= theta.colwise().mean();
  const Eigen::RowVectorXd proportions = y.colwise().mean();
  const double reference_proportion = 1.0 - proportions.sum();
  double step = 1.0;

  for (Eigen::Index i = 0; i < theta.rows(); ++i) {
    double target_reference = 1.0;
    for (Eigen::Index k = 0; k < theta.cols(); ++k) {
      const double target = y(i, k) + theta(i, k);
      const double direction = target - proportions(k);
      if (direction < 0.0) {
        step = std::min(step, proportions(k) / -direction);
      }
      target_reference -= target;
    }
    const double reference_direction = target_reference - reference_proportion;
    if (reference_direction < 0.0) {
      step = std::min(step, reference_proportion / -reference_direction);
    }
  }

  step = std::clamp(step, 0.0, 1.0);
  if (step > 0.0 && step < 1.0) {
    step *= 1.0 - std::sqrt(std::numeric_limits<double>::epsilon());
  }

  for (Eigen::Index i = 0; i < theta.rows(); ++i) {
    for (Eigen::Index k = 0; k < theta.cols(); ++k) {
      const double anchor = proportions(k) - y(i, k);
      theta(i, k) = anchor + step * (theta(i, k) - anchor);
    }
  }

  return theta;
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
