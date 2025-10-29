/**
 * @file
 * @brief Scoring metrics for model evaluation
 */

#pragma once

#include "losses/loss.h"
#include <Eigen/Core>
#include <memory>

namespace slope {

/**
 * Computes the Area Under the ROC Curve (AUC-ROC) for binary classification.
 *
 * @param scores Vector of prediction scores/probabilities for each sample
 * @param labels Vector of true binary labels (0 or 1)
 * @return The AUC-ROC score between 0 and 1, where 1 indicates perfect
 * prediction and 0.5 indicates random prediction
 *
 * Both input vectors must have the same length. The function calculates how
 * well the scores discriminate between the positive and negative classes.
 */
double
binaryRocAuc(const Eigen::VectorXd& scores, const Eigen::VectorXd& labels);

/**
 * Computes the Area Under the ROC Curve (AUC-ROC) for binary and multi-class
 * classification.
 *
 * @param scores Matrix of prediction scores/probabilities (samples x classes)
 * @param labels Matrix of true labels in one-hot encoded format (samples x
 * classes)
 * @return The average AUC-ROC score across all classes, between 0 and 1
 *
 * For binary classification, this reduces to standard binary AUC-ROC.
 * For multi-class, computes one-vs-rest AUC-ROC for each class and averages.
 * Both input matrices must have the same dimensions.
 */
double
rocAuc(const Eigen::MatrixXd& scores, const Eigen::MatrixXd& labels);

/**
 * @brief Base class for scoring metrics used in regularized generalized linear
 * regression.
 *
 * This abstract class defines the interface for computing various performance
 * metrics that evaluate model predictions against true responses. Derived
 * classes implement specific metrics like MSE, MAE, AUC-ROC etc.
 */
class Score
{
public:
  virtual ~Score() = default;

  /**
   * Determines if a score value indicates worse performance.
   * @param other The score to compare against
   * @param current The current score
   * @return true if other score is worse than current
   */
  virtual bool isWorse(double other, double current) const = 0;

  /**
   * Returns the initial/default value for this scoring metric.
   * @return Initial score value
   */
  virtual double initValue() const = 0;

  /**
   * Returns a comparator function for comparing score values.
   * @return Function object that compares two score values
   */
  std::function<bool(double, double)> getComparator() const;

  /**
   * Evaluates the scoring metric given predictions and true responses.
   * @param eta Matrix of model predictions
   * @param y Matrix of true responses
   * @param loss Loss function used in the model
   * @return Computed score value
   */
  virtual double eval(const Eigen::MatrixXd& eta,
                      const Eigen::MatrixXd& y,
                      const std::unique_ptr<Loss>& loss) const = 0;

  /**
   * Factory method to create specific Score implementations.
   * @param metric Name of the scoring metric to create
   * @return Unique pointer to created Score object
   */
  static std::unique_ptr<Score> create(const std::string& metric);
};

/**
 * @brief Scoring metric that aims to minimize the score value.
 *
 * This class implements a scoring metric where lower values indicate better
 * performance (e.g., MSE, MAE, deviance). It overrides the base Score class
 * comparison methods to establish that higher values are worse.
 */
class MinimizeScore : public Score
{
public:
  /**
   * Determines if score a is worse than score b.
   * @param a First score to compare
   * @param b Second score to compare
   * @return true if a > b (higher values are worse)
   */
  bool isWorse(double a, double b) const override;

  /**
   * Returns the initial value for minimization metrics.
   * @return Positive infinity as the worst possible score
   */
  double initValue() const override;
};

/**
 * @brief Scoring metric that aims to maximize the score value.
 *
 * This class implements a scoring metric where higher values indicate better
 * performance (e.g., R², AUC-ROC, accuracy). It overrides the base Score class
 * comparison methods to establish that lower values are worse.
 */
class MaximizeScore : public Score
{
public:
  /**
   * Determines if score a is worse than score b.
   * @param a First score to compare
   * @param b Second score to compare
   * @return true if a < b (lower values are worse)
   */
  bool isWorse(double a, double b) const override;

  /**
   * Returns the initial value for maximization metrics.
   * @return Negative infinity as the worst possible score
   */
  double initValue() const override;
};

/**
 * @brief Mean Squared Error (MSE) scoring metric.
 *
 * Computes the average squared difference between predictions and true
 * responses. Inherits from MinimizeScore since lower MSE values indicate better
 * fit.
 *
 * MSE = (1/n) Σ(y_i - η_i)²
 * where:
 * - n is the number of samples
 * - y_i is the true response for sample i
 * - η_i is the predicted response for sample i
 */
class MSE : public MinimizeScore
{
public:
  /**
   * Evaluates the MSE between predictions and true responses.
   * @param eta Matrix of model predictions
   * @param y Matrix of true responses
   * @return Computed MSE value
   */
  double eval(const Eigen::MatrixXd& eta,
              const Eigen::MatrixXd& y,
              const std::unique_ptr<Loss>&) const override;
};

/**
 * @brief Mean Absolute Error (MAE) scoring metric.
 *
 * Computes the average absolute difference between predictions and true
 * responses. Inherits from MinimizeScore since lower MAE values indicate better
 * fit.
 *
 * MAE = (1/n) Σ|y_i - η_i|
 * where:
 * - n is the number of samples
 * - y_i is the true response for sample i
 * - η_i is the predicted response for sample i
 *
 * MAE is more robust to outliers compared to MSE as it uses absolute
 * rather than squared differences.
 */
class MAE : public MinimizeScore
{
public:
  /**
   * Evaluates the MAE between predictions and true responses.
   * @param eta Matrix of model predictions
   * @param y Matrix of true responses
   * @return Computed MAE value
   */
  double eval(const Eigen::MatrixXd& eta,
              const Eigen::MatrixXd& y,
              const std::unique_ptr<Loss>&) const override;
};

/**
 * @brief Classification Accuracy scoring metric.
 *
 * Computes the proportion of correct predictions in classification problems.
 * Inherits from MaximizeScore since higher accuracy values indicate better
 * performance.
 *
 * Accuracy = (number of correct predictions) / (total number of predictions)
 *
 * For binary classification, predictions are thresholded at 0.5.
 * For multi-class problems, the class with highest probability is selected
 * as the prediction for each sample.
 */
class Accuracy : public MaximizeScore
{
public:
  /**
   * Evaluates the classification accuracy between predictions and true labels.
   * @param eta Matrix of model predictions/probabilities
   * @param y Matrix of true labels (one-hot encoded)
   * @param loss Loss function used to transform predictions if needed
   * @return Proportion of correct predictions (between 0 and 1)
   */
  double eval(const Eigen::MatrixXd& eta,
              const Eigen::MatrixXd& y,
              const std::unique_ptr<Loss>& loss) const override;
};

/**
 * @brief Misclassification Rate scoring metric.
 *
 * Computes the proportion of incorrect predictions in classification problems.
 * Inherits from MinimizeScore since lower misclassification rates indicate
 * better performance.
 *
 * Misclassification Rate = 1 - Accuracy
 *                       = (number of incorrect predictions) / (total number of
 * predictions)
 *
 * For binary classification, predictions are thresholded at 0.5.
 * For multi-class problems, the class with highest probability is selected
 * as the prediction for each sample.
 */
class MisClass : public MinimizeScore
{
public:
  /**
   * Evaluates the misclassification rate between predictions and true labels.
   * @param eta Matrix of model predictions/probabilities
   * @param y Matrix of true labels (one-hot encoded)
   * @param loss Loss function used to transform predictions if needed
   * @return Proportion of incorrect predictions (between 0 and 1)
   */
  double eval(const Eigen::MatrixXd& eta,
              const Eigen::MatrixXd& y,
              const std::unique_ptr<Loss>& loss) const override;
};

/**
 * @brief Deviance scoring metric.
 *
 * Computes the statistical deviance, which is -2 times the log-likelihood
 * of the predictions under the specified loss function. Inherits from
 * MinimizeScore since lower deviance indicates better model fit.
 *
 * Deviance = 2 * ( log L(saturated) - log L (fitted) )
 *
 * The exact form of deviance depends on the loss function used.
 */
class Deviance : public MinimizeScore
{
public:
  /**
   * Evaluates the deviance between predictions and true responses.
   * @param eta Matrix of model predictions
   * @param y Matrix of true responses
   * @param loss Loss function that defines the likelihood
   * @return Computed deviance value
   */
  double eval(const Eigen::MatrixXd& eta,
              const Eigen::MatrixXd& y,
              const std::unique_ptr<Loss>& loss) const override;
};

/**
 * @brief Area Under the ROC Curve (AUC-ROC) scoring metric.
 *
 * Computes the area under the Receiver Operating Characteristic curve,
 * which plots the true positive rate against false positive rate at
 * various classification thresholds. Inherits from MaximizeScore since
 * higher AUC values indicate better discriminative ability.
 *
 * For binary classification:
 * - AUC of 1.0 represents perfect prediction
 * - AUC of 0.5 represents random prediction
 * - AUC of 0.0 represents perfectly inverse prediction
 *
 * For multi-class problems:
 * - Computes the one-vs-rest AUC for each class
 * - Returns the average across all classes
 */
class AUC : public MaximizeScore
{
public:
  /**
   * Evaluates the AUC-ROC score for predictions.
   * @param eta Matrix of model predictions/probabilities
   * @param y Matrix of true labels (one-hot encoded)
   * @param loss Loss function used to transform predictions if needed
   * @return AUC-ROC score (between 0 and 1)
   */
  double eval(const Eigen::MatrixXd& eta,
              const Eigen::MatrixXd& y,
              const std::unique_ptr<Loss>& loss) const override;
};

} // namespace slope
