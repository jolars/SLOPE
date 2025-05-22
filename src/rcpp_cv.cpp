#include "setup_model.h"
#include "slope/cv.h"
#include "slope/slope.h"
#include <RcppEigen.h>

namespace Rcpp {
template<>
SEXP
wrap(const slope::GridResult& gr)
{
  Rcpp::List result =
    Rcpp::List::create(Rcpp::Named("score") = Rcpp::wrap(gr.score),
                       Rcpp::Named("params") = Rcpp::wrap(gr.params),
                       Rcpp::Named("alphas") = Rcpp::wrap(gr.alphas),
                       Rcpp::Named("mean_scores") = Rcpp::wrap(gr.mean_scores),
                       Rcpp::Named("std_errors") = Rcpp::wrap(gr.std_errors));
  return result;
}
}

template<typename T>
Rcpp::List
cvImpl(T& x,
       const Eigen::MatrixXd& y,
       const Rcpp::List& cv_args,
       const Rcpp::List& model_args)
{
  using Eigen::ArrayXd;
  using Eigen::VectorXd;
  using Rcpp::as;
  using Rcpp::Named;
  using Rcpp::wrap;

  slope::Slope model = setupModel(model_args);

  auto lambda_type = as<std::string>(model_args["lambda_type"]);

  if (lambda_type != "user") {
    model.setLambdaType(lambda_type);
  }

  auto cv_config = slope::CvConfig();

  std::map<std::string, std::vector<double>> hyperparams;
  hyperparams["q"] = as<std::vector<double>>(cv_args["q"]);
  hyperparams["gamma"] = as<std::vector<double>>(cv_args["gamma"]);

  cv_config.hyperparams = hyperparams;
  cv_config.metric = as<std::string>(cv_args["metric"]);
  cv_config.predefined_folds =
    as<std::vector<std::vector<std::vector<int>>>>(cv_args["predefined_folds"]);

  auto res = crossValidate(model, x, y, cv_config);

  return Rcpp::List::create(Named("results") = wrap(res.results),
                            Named("best_params") = wrap(res.best_params),
                            Named("best_score") = wrap(res.best_score),
                            Named("best_ind") = wrap(res.best_ind),
                            Named("best_alpha_ind") = wrap(res.best_alpha_ind));
}

// [[Rcpp::export]]
Rcpp::List
cvSparseCpp(Eigen::SparseMatrix<double>& x,
            const Eigen::MatrixXd& y,
            const Rcpp::List& cv_args,
            const Rcpp::List& model_args)
{
  return cvImpl(x, y, cv_args, model_args);
}

// [[Rcpp::export]]
Rcpp::List
cvDenseCpp(Eigen::MatrixXd& x,
           const Eigen::MatrixXd& y,
           const Rcpp::List& cv_args,
           const Rcpp::List& model_args)
{
  if (!x.allFinite()) {
    throw std::invalid_argument("Input matrix must not contain missing values");
  }
  return cvImpl(x, y, cv_args, model_args);
}
