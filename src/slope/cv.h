/**
 * @file
 * @brief Cross-validation functionality for SLOPE models
 *
 * This file provides data structures and functions for performing k-fold
 * cross-validation on SLOPE models to find optimal hyperparameters.
 * It includes functionality for parameter grid searches, fold generation,
 * and result aggregation with support for parallel processing.
 */

#pragma once

#include "folds.h"
#include "score.h"
#include "slope.h"
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace slope {

/**
 * @brief Stores cross-validation results for a specific set of hyperparameters
 *
 * This struct contains evaluation scores, parameters, and statistics for a
 * single hyperparameter configuration across all cross-validation folds.
 */
struct GridResult
{
  /// Matrix of evaluation scores indexed by (fold, alpha) where each row
  /// represents a fold and each column represents an alpha value
  /// (regularization weight)
  Eigen::MatrixXd score;

  /// Map of hyperparameter names to their values for the configuration
  std::map<std::string, double> params;

  /// Array of regularization parameters used in the regularization path
  Eigen::ArrayXd alphas;

  /// Array of scores averaged across all folds for each alpha value
  Eigen::ArrayXd mean_scores;

  /// Array of standard errors of the scores across folds for each alpha value,
  /// useful for estimating score variability
  Eigen::ArrayXd std_errors;
};

/**
 * @brief Contains overall results from a cross-validation process
 *
 * This struct aggregates results from cross-validation across multiple
 * hyperparameter combinations, including information about the optimal
 * configuration.
 */
struct CvResult
{
  /// Vector of GridResult objects containing performance metrics for each
  /// hyperparameter configuration evaluated
  std::vector<GridResult> results;

  /// Map of hyperparameter names to their optimal values based on the
  /// cross-validation results
  std::map<std::string, double> best_params;

  /// The score achieved by the optimal hyperparameter configuration
  double best_score;

  /// Index of the best performing configuration in the results vector
  int best_ind;

  /// Index of the optimal alpha value within the regularization path for the
  /// best configuration
  int best_alpha_ind;
};

/**
 * @brief Configuration settings for cross-validation
 *
 * This struct specifies the parameters used to control the cross-validation
 * process, including fold count, evaluation metric, random seed, and
 * hyperparameter grid.
 */
struct CvConfig
{
  /// Number of folds for cross-validation (default: 10)
  int n_folds = 10;

  /// Number of times to repeat the cross-validation (default: 1)
  int n_repeats = 1;

  /// Evaluation metric used for model assessment (default: "mse")
  std::string metric = "mse";

  /// Seed for random number generator to ensure reproducibility (default: 42)
  uint64_t random_seed = 42;

  /// Map of hyperparameter names to vectors of values to evaluate
  std::map<std::string, std::vector<double>> hyperparams;

  /// Map of hyperparameter names to vectors of values to evaluate
  std::map<std::string, std::vector<double>> default_hyperparams = {
    { "q", { 0.1 } },
    { "gamma", { 0.0 } },
  };

  /// Optional user-defined fold assignments for custom cross-validation splits
  std::optional<std::vector<std::vector<std::vector<int>>>> predefined_folds;
};

/**
 * @brief Creates a grid of parameter combinations from parameter value ranges
 *
 * @param param_values Map of parameter names to vectors of possible values
 * @return std::vector<std::map<std::string, double>> Vector of parameter
 * combinations
 *
 * This function takes a map where keys are parameter names and values are
 * vectors of possible values for each parameter, then generates all possible
 * combinations.
 */
std::vector<std::map<std::string, double>>
createGrid(const std::map<std::string, std::vector<double>>& param_values);

/**
 * @brief Identifies the best parameters from cross-validation results
 *
 * @param cv_result Cross-validation results to analyze and update
 * @param scorer The scoring metric used to evaluate performance
 *
 * This function examines all cross-validation results across the parameter grid
 * and updates the cv_result with information about the best parameter
 * combination, including the best score and corresponding indices.
 */
void
findBestParameters(CvResult& cv_result, const std::unique_ptr<Score>& scorer);

/**
 * @brief Performs cross-validation on a SLOPE model to select optimal
 * hyperparameters
 *
 * @tparam MatrixType Type of design matrix (supports both dense and sparse
 * matrices)
 * @param model The SLOPE model to be cross-validated
 * @param x The design matrix containing predictors
 * @param y_in The response matrix
 * @param config Configuration parameters for cross-validation (optional)
 * @return CvResult Object containing cross-validation results, including best
 * parameters, regularization paths, and performance metrics across folds
 *
 * This function implements k-fold cross-validation for SLOPE models, evaluating
 * model performance across a grid of hyperparameters. For each hyperparameter
 * combination, the function:
 * 1. Splits the data into training and validation sets according to the fold
 * configuration
 * 2. Fits the model on each training set and evaluates on the validation set
 * 3. Computes the specified evaluation metric for each regularization parameter
 * 4. Averages results across folds to select optimal hyperparameters
 *
 * The function supports parallel processing with OpenMP when available.
 */
template<typename MatrixType>
CvResult
crossValidate(Slope model,
              MatrixType& x,
              const Eigen::MatrixXd& y_in,
              const CvConfig& config = CvConfig())
{
  CvResult cv_result;

  int n = y_in.rows();

  auto loss = setupLoss(model.getLossType());

  auto y = loss->preprocessResponse(y_in);
  auto scorer = Score::create(config.metric);

  auto hyperparams = config.default_hyperparams;

  // Override with user-specified parameters
  for (const auto& [key, values] : config.hyperparams) {
    hyperparams[key] = values;
  }

  auto grid = createGrid(hyperparams);

  // Total number of evaluations (n_repeats * n_folds)
  Folds folds =
    config.predefined_folds.has_value()
      ? Folds(config.predefined_folds.value())
      : Folds(n, config.n_folds, config.n_repeats, config.random_seed);

  int n_evals = folds.numEvals();

  for (const auto& params : grid) {
    GridResult result;
    result.params = params;

    double q = params.at("q");
    double gamma = params.at("gamma");

    model.setQ(q);

    auto initial_path = model.path(x, y);

    result.alphas = initial_path.getAlpha();
    int n_alpha = result.alphas.size();

    assert((result.alphas > 0).all());

    Eigen::MatrixXd scores = Eigen::MatrixXd::Zero(n_evals, n_alpha);

    Eigen::setNbThreads(1);

    // Thread-safety for exceptions
    std::vector<std::string> thread_errors(n_evals);
    bool had_exception = false;

#ifdef _OPENMP
    omp_set_max_active_levels(1);
#pragma omp parallel for num_threads(Threads::get())                           \
  shared(scores, thread_errors, had_exception)
#endif
    for (int i = 0; i < n_evals; ++i) {
      try {
        auto [rep, fold] = std::div(i, folds.numFolds());

        Slope thread_model = model;
        thread_model.setModifyX(true);

        // TODO: Maybe consider not copying at all?
        auto [x_train, y_train, x_test, y_test] = folds.split(x, y, fold, rep);

        auto path = thread_model.path(x_train, y_train, result.alphas);

        if (gamma > 0) {
          path = thread_model.relax(path, x_train, y_train, gamma);
        }

        for (int j = 0; j < n_alpha; ++j) {
          auto eta = path(j).predict(x_test, "linear");
          scores(i, j) = scorer->eval(eta, y_test, loss);
        }
      } catch (const std::exception& e) {
        thread_errors[i] = e.what();
#ifdef _OPENMP
#pragma omp atomic write
#endif
        had_exception = true;
      } catch (...) {
        thread_errors[i] = "Unknown exception";
#ifdef _OPENMP
#pragma omp atomic write
#endif
        had_exception = true;
      }
    }

    if (had_exception) {
      std::string error_message = "Exception(s) during cross-validation:\n";
      for (int i = 0; i < n_evals; ++i) {
        if (!thread_errors[i].empty()) {
          error_message +=
            "Fold " + std::to_string(i) + ": " + thread_errors[i] + "\n";
        }
      }
      throw std::runtime_error(error_message);
    }

    result.mean_scores = scores.colwise().mean();
    result.std_errors = stdDevs(scores).array() / std::sqrt(n_evals);
    result.score = std::move(scores);
    cv_result.results.push_back(result);
  }

#ifdef _OPENMP
  Eigen::setNbThreads(0);
#endif

  findBestParameters(cv_result, scorer);

  return cv_result;
}

} // namespace slope
