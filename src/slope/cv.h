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
  Eigen::MatrixXd score;                // indexed by (fold, alpha)
  std::map<std::string, double> params; // hyperparams (q, etc.)
  Eigen::ArrayXd alphas;                // the sequence of alphas from the path
  Eigen::ArrayXd
    mean_scores; // averaged over folds for each (param,alpha) combo
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
  std::vector<GridResult> results;

  std::map<std::string, double> best_params;
  double best_score;
  int best_ind;
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
  int n_folds = 10;
  std::string metric = "mse";
  uint64_t random_seed = 42;
  std::map<std::string, std::vector<double>> hyperparams = { { "q", { 0.1 } } };
  std::optional<std::vector<std::vector<int>>> predefined_folds;
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
  auto grid = createGrid(config.hyperparams);

  Folds folds = config.predefined_folds.has_value()
                  ? Folds(config.predefined_folds.value())
                  : Folds(n, config.n_folds, config.random_seed);

  for (const auto& params : grid) {
    GridResult result;
    result.params = params;
    model.setQ(params.at("q"));

    auto initial_path = model.path(x, y);
    result.alphas = initial_path.getAlpha();

    assert((result.alphas > 0).all());

    Eigen::MatrixXd scores =
      Eigen::MatrixXd::Zero(config.n_folds, result.alphas.size());

    Eigen::setNbThreads(1);

#ifdef _OPENMP
    omp_set_max_active_levels(1);
#pragma omp parallel for num_threads(Threads::get()) shared(scores)
#endif
    for (int i = 0; i < config.n_folds; ++i) {
      Slope thread_model = model;
      thread_model.setModifyX(true);

      // TODO: Maybe consider not copying at all?
      auto [x_train, y_train, x_test, y_test] = folds.split(x, y, i);
      auto path = thread_model.path(x_train, y_train, result.alphas);

      for (int j = 0; j < result.alphas.size(); ++j) {
        auto eta = path(j).predict(x_test, "linear");
        scores(i, j) = scorer->eval(eta, y_test, loss);
      }
    }

    result.mean_scores = scores.colwise().mean();
    result.std_errors = stdDevs(scores).array() / std::sqrt(config.n_folds);
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
