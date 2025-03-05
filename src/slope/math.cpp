#include "math.h"
#include "constants.h"

namespace slope {

Eigen::VectorXd
logSumExp(const Eigen::MatrixXd& a)
{
  Eigen::VectorXd max_vals = a.rowwise().maxCoeff();
  Eigen::ArrayXd sum_exp =
    (a.colwise() - max_vals).array().exp().rowwise().sum();

  return max_vals.array() + sum_exp.max(constants::P_MIN).log();
}

Eigen::MatrixXd
softmax(const Eigen::MatrixXd& a)
{
  Eigen::VectorXd shift = a.rowwise().maxCoeff();
  Eigen::MatrixXd exp_a = (a.colwise() - shift).array().exp();
  Eigen::ArrayXd row_sums = exp_a.rowwise().sum();
  return exp_a.array().colwise() / row_sums;
}

std::vector<int>
setUnion(const std::vector<int>& a, const std::vector<int>& b)
{
  std::vector<int> out;
  std::set_union(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  return out;
}

std::vector<int>
setDiff(const std::vector<int>& a, const std::vector<int>& b)
{
  std::vector<int> out;
  std::set_difference(
    a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));

  return out;
}

Eigen::ArrayXd
geomSpace(const double start, const double end, const int n)
{
  if (n == 1) {
    return Eigen::ArrayXd::Constant(1, start);
  } else {
    return Eigen::ArrayXd::LinSpaced(n, std::log(start), std::log(end)).exp();
  }
}

Eigen::VectorXd
l2Norms(const Eigen::SparseMatrix<double>& x)
{
  const int p = x.cols();

  Eigen::VectorXd out(p);

  for (int j = 0; j < p; ++j) {
    out(j) = x.col(j).norm();
  }

  return out;
}

Eigen::VectorXd
l2Norms(const Eigen::MatrixXd& x)
{
  return x.colwise().norm();
}

Eigen::VectorXd
means(const Eigen::SparseMatrix<double>& x)
{
  const int n = x.rows();
  const int p = x.cols();

  Eigen::VectorXd out(p);

  for (int j = 0; j < p; ++j) {
    out(j) = x.col(j).sum() / n;
  }

  return out;
}

Eigen::VectorXd
means(const Eigen::MatrixXd& x)
{
  return x.colwise().mean();
}

Eigen::VectorXd
ranges(const Eigen::SparseMatrix<double>& x)
{
  const int p = x.cols();

  Eigen::VectorXd out(p);

  for (int j = 0; j < p; ++j) {
    double x_j_max = 0.0;
    double x_j_min = 0.0;

    for (typename Eigen::SparseMatrix<double>::InnerIterator it(x, j); it;
         ++it) {
      x_j_max = std::max(x_j_max, it.value());
      x_j_min = std::min(x_j_min, it.value());
    }

    out(j) = x_j_max - x_j_min;
  }

  return out;
}

Eigen::VectorXd
ranges(const Eigen::MatrixXd& x)
{
  return x.colwise().maxCoeff() - x.colwise().minCoeff();
}

Eigen::VectorXd
maxAbs(const Eigen::SparseMatrix<double>& x)
{
  const int p = x.cols();

  Eigen::VectorXd out(p);

  for (int j = 0; j < p; ++j) {
    double x_j_maxabs = 0.0;

    for (typename Eigen::SparseMatrix<double>::InnerIterator it(x, j); it;
         ++it) {
      x_j_maxabs = std::max(x_j_maxabs, std::abs(it.value()));
    }

    out(j) = x_j_maxabs;
  }

  return out;
}

Eigen::VectorXd
maxAbs(const Eigen::MatrixXd& x)
{
  return x.cwiseAbs().colwise().maxCoeff();
}

Eigen::VectorXd
mins(const Eigen::SparseMatrix<double>& x)
{
  const int p = x.cols();

  Eigen::VectorXd out(p);

  for (int j = 0; j < p; ++j) {
    double x_j_min = 0.0;

    for (typename Eigen::SparseMatrix<double>::InnerIterator it(x, j); it;
         ++it) {
      x_j_min = std::min(x_j_min, it.value());
    }

    out(j) = x_j_min;
  }

  return out;
}

Eigen::VectorXd
mins(const Eigen::MatrixXd& x)
{
  return x.colwise().minCoeff();
}

} // namespace slope
