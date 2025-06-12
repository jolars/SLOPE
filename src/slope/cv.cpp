#include "cv.h"

namespace slope {

std::vector<std::map<std::string, double>>
createGrid(const std::map<std::string, std::vector<double>>& param_values)
{
  std::vector<std::map<std::string, double>> grid;

  if (param_values.empty()) {
    return grid;
  }

  // Start with first parameter
  auto it = param_values.begin();
  for (double value : it->second) {
    std::map<std::string, double> point;
    point[it->first] = value;
    grid.push_back(point);
  }

  // Add remaining parameters
  for (++it; it != param_values.end(); ++it) {
    std::vector<std::map<std::string, double>> new_grid;
    for (const auto& existing_point : grid) {
      for (double value : it->second) {
        auto new_point = existing_point;
        new_point[it->first] = value;
        new_grid.push_back(new_point);
      }
    }
    grid = std::move(new_grid);
  }

  return grid;
}

void
findBestParameters(CvResult& cv_result, const std::unique_ptr<Score>& scorer)
{
  cv_result.best_score = scorer->initValue();
  auto comp = scorer->getComparator();

  for (size_t i = 0; i < cv_result.results.size(); ++i) {
    auto result = cv_result.results[i];
    int best_alpha_ind = whichBest(result.mean_scores, comp);
    double current_score = result.mean_scores(best_alpha_ind);

    assert(best_alpha_ind >= 0 && best_alpha_ind < result.alphas.size());
    assert(result.alphas(best_alpha_ind) > 0);

    if (scorer->isWorse(cv_result.best_score, current_score)) {
      cv_result.best_ind = i;
      cv_result.best_score = current_score;
      cv_result.best_params = result.params;
      cv_result.best_alpha_ind = best_alpha_ind;

      cv_result.best_params["alpha"] = result.alphas(best_alpha_ind);
    }
  }
}

} // namespace slope
