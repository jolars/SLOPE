#include "score.h"
#include "constants.h"
#include <Eigen/Core>
#include <memory>
#include <vector>

namespace slope {

double
binaryRocAuc(const Eigen::VectorXd& scores, const Eigen::VectorXd& labels)
{
  int n = scores.rows();
  std::vector<std::pair<double, bool>> pairs(n);

  for (int i = 0; i < n; ++i) {
    pairs[i] = { scores(i), labels(i) > 0.5 };
  }

  std::sort(pairs.begin(), pairs.end(), [](const auto& a, const auto& b) {
    return a.first > b.first;
  });

  int pos_count = labels.sum();
  int neg_count = n - pos_count;

  if (pos_count == 0 || neg_count == 0)
    return 0.5;

  double auc = 0.0;
  double fp = 0, tp = 0;
  double fp_prev = 0, tp_prev = 0;

  for (size_t i = 0; i < pairs.size(); ++i) {
    if (pairs[i].second) {
      tp++;
    } else {
      fp++;
    }

    if (i < pairs.size() - 1 && pairs[i].first != pairs[i + 1].first) {
      double fp_rate = fp / neg_count;
      double tp_rate = tp / pos_count;
      auc += (fp_rate - fp_prev) * (tp_rate + tp_prev) * 0.5;
      fp_prev = fp_rate;
      tp_prev = tp_rate;
    }
  }

  auc += (1 - fp_prev) * (1 + tp_prev) * 0.5;
  return auc;
}

double
rocAuc(const Eigen::MatrixXd& scores, const Eigen::MatrixXd& labels)
{
  int n = scores.rows();
  int m = scores.cols(); // number of classes

  // Validate dimensions
  if (scores.rows() != labels.rows() || scores.cols() != labels.cols()) {
    throw std::invalid_argument("scores and labels dimensions must match");
  }

  if (m == 1) {
    // Binary classification case
    return binaryRocAuc(scores, labels);
  }

  // Compute macro-average AUC (one-vs-rest)
  double total_auc = 0.0;
  int valid_classes = 0;

  for (int c = 0; c < m; ++c) {
    // Extract scores for current class
    Eigen::VectorXd class_scores = scores.col(c);

    // Create binary labels for current class (1 for current class, 0 for
    // others)
    Eigen::VectorXd binary_labels = Eigen::VectorXd::Zero(n);
    for (int i = 0; i < n; ++i) {
      binary_labels(i) = (labels(i, c) > 0.5) ? 1.0 : 0.0;
    }

    // Skip if class has no positive or negative samples
    if (binary_labels.sum() == 0 || binary_labels.sum() == n) {
      continue;
    }

    total_auc += binaryRocAuc(class_scores, binary_labels);
    valid_classes++;
  }

  return valid_classes > 0 ? total_auc / valid_classes : 0.5;
}

std::function<bool(double, double)>
Score::getComparator() const
{
  return [this](double a, double b) { return this->isWorse(a, b); };
}

bool
MinimizeScore::isWorse(double a, double b) const
{
  return a > b;
}

double
MinimizeScore::initValue() const
{
  return constants::POS_INF;
}

bool
MaximizeScore::isWorse(double a, double b) const
{
  return a < b;
}

double
MaximizeScore::initValue() const
{
  return constants::NEG_INF;
}

double
MSE::eval(const Eigen::MatrixXd& eta,
          const Eigen::MatrixXd& y,
          const std::unique_ptr<Loss>&) const
{
  return (y - eta).squaredNorm() / y.rows();
}

double
MAE::eval(const Eigen::MatrixXd& eta,
          const Eigen::MatrixXd& y,
          const std::unique_ptr<Loss>&) const
{
  return (y - eta).cwiseAbs().mean();
}

double
Accuracy::eval(const Eigen::MatrixXd& eta,
               const Eigen::MatrixXd& y,
               const std::unique_ptr<Loss>& loss) const
{
  Eigen::MatrixXd y_pred = loss->predict(eta);

  // There is a bug in Eigen, so we need to cast to double first.
  Eigen::MatrixXd comparison = (y_pred.array() == y.array()).cast<double>();
  return comparison.sum() / y.rows();
}

double
MisClass::eval(const Eigen::MatrixXd& eta,
               const Eigen::MatrixXd& y,
               const std::unique_ptr<Loss>& loss) const
{
  Eigen::MatrixXd y_pred = loss->predict(eta);
  Eigen::MatrixXd comparison = (y_pred.array() == y.array()).cast<double>();
  return 1 - comparison.sum() / y.rows();
}

double
Deviance::eval(const Eigen::MatrixXd& eta,
               const Eigen::MatrixXd& y,
               const std::unique_ptr<Loss>& loss) const
{
  return loss->deviance(eta, y);
}

double
AUC::eval(const Eigen::MatrixXd& eta,
          const Eigen::MatrixXd& y,
          const std::unique_ptr<Loss>& loss) const
{
  Eigen::MatrixXd probs = loss->inverseLink(eta);
  return rocAuc(probs, y);
}

std::unique_ptr<Score>
Score::create(const std::string& metric)
{
  if (metric == "mse")
    return std::make_unique<MSE>();
  if (metric == "mae")
    return std::make_unique<MAE>();
  if (metric == "accuracy")
    return std::make_unique<Accuracy>();
  if (metric == "misclass")
    return std::make_unique<MisClass>();
  if (metric == "deviance")
    return std::make_unique<Deviance>();
  if (metric == "auc")
    return std::make_unique<AUC>();

  throw std::invalid_argument("Unknown metric: " + metric);
}

};
