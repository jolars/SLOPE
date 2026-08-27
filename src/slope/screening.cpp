/**
 * @file
 * @brief Implementation of screening rules for SLOPE regression optimization
 */

#include "kkt_check.h"
#include <Eigen/Core>
#include <cassert>
#include <numeric>
#include <slope/math.h>
#include <slope/screening.h>
#include <slope/utils.h>
#include <stdexcept>

namespace slope {

typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;
typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXXb;

std::vector<int>
activeSet(const Eigen::VectorXd& beta)
{
  ArrayXb active = beta.array() != 0.0;

  return which(active);
}

std::vector<int>
strongSet(const Eigen::VectorXd& gradient_prev,
          const Eigen::ArrayXd& lambda,
          const Eigen::ArrayXd& lambda_prev)
{
  using Eigen::VectorXd;
  using Eigen::VectorXi;

  int pm = gradient_prev.size();

  assert(lambda_prev.size() == lambda.size() &&
         "lambda_prev and lambda must have the same length");
  assert((lambda <= lambda_prev).all() &&
         "New lambda values must be smaller than or equal to previous values");

  const VectorXd abs_grad = gradient_prev.reshaped().cwiseAbs();
  std::vector<int> ord = sortIndex(abs_grad, true);

  assert(abs_grad.size() == lambda.size());

  const VectorXd tmp =
    abs_grad(ord).array().eval() + lambda_prev - 2.0 * lambda;

  int i = 0;
  int k = 0;

  double s = 0;

  while (i + k < pm) {
    s += tmp(k + i);

    if (s >= 0) {
      k = k + i + 1;
      i = 0;
      s = 0;
    } else {
      i++;
    }
  }

  ArrayXb active_set = ArrayXb::Zero(pm);
  active_set.head(k).setOnes();

  // restore order
  inversePermute(active_set, ord);

  return which(active_set);
}

// NoScreening implementation
std::vector<int>
NoScreening::initialize(int feature_count, int)
{
  std::vector<int> working_set(feature_count);
  std::iota(working_set.begin(), working_set.end(), 0);
  return working_set;
}

void
NoScreening::screen(std::vector<int>&,
                    Eigen::VectorXd&,
                    const Eigen::ArrayXd&,
                    const Eigen::ArrayXd&,
                    const Eigen::VectorXd&)
{
}

bool
NoScreening::checkKktViolations(Eigen::VectorXd&,
                                const Eigen::VectorXd&,
                                const Eigen::ArrayXd&,
                                std::vector<int>&,
                                const Eigen::MatrixXd&,
                                const Eigen::MatrixXd&,
                                const Eigen::VectorXd&,
                                const Eigen::VectorXd&,
                                JitNormalization)
{
  return true;
}

bool
NoScreening::checkKktViolations(Eigen::VectorXd&,
                                const Eigen::VectorXd&,
                                const Eigen::ArrayXd&,
                                std::vector<int>&,
                                const Eigen::SparseMatrix<double>&,
                                const Eigen::MatrixXd&,
                                const Eigen::VectorXd&,
                                const Eigen::VectorXd&,
                                JitNormalization)
{
  return true;
}

bool
NoScreening::checkKktViolations(Eigen::VectorXd&,
                                const Eigen::VectorXd&,
                                const Eigen::ArrayXd&,
                                std::vector<int>&,
                                const Eigen::Map<Eigen::MatrixXd>&,
                                const Eigen::MatrixXd&,
                                const Eigen::VectorXd&,
                                const Eigen::VectorXd&,
                                JitNormalization)
{
  return true;
}

bool
NoScreening::checkKktViolations(Eigen::VectorXd&,
                                const Eigen::VectorXd&,
                                const Eigen::ArrayXd&,
                                std::vector<int>&,
                                const Eigen::Map<Eigen::SparseMatrix<double>>&,
                                const Eigen::MatrixXd&,
                                const Eigen::VectorXd&,
                                const Eigen::VectorXd&,
                                JitNormalization)
{
  return true;
}

std::string
NoScreening::toString() const
{
  return "none";
}

// StrongScreening implementation
std::vector<int>
StrongScreening::initialize(int, int alpha_max_ind)
{
  return { alpha_max_ind };
}

void
StrongScreening::screen(std::vector<int>& working_set,
                        Eigen::VectorXd& gradient,
                        const Eigen::ArrayXd& lambda_curr,
                        const Eigen::ArrayXd& lambda_prev,
                        const Eigen::VectorXd& beta)
{
  if (lambda_curr(0) == 0.0) {
    working_set.resize(beta.size());
    std::iota(working_set.begin(), working_set.end(), 0);
    return;
  }

  std::vector<int> active_set = activeSet(beta);
  strong_set = strongSet(gradient, lambda_curr, lambda_prev);
  strong_set = setUnion(strong_set, active_set);

  working_set = setUnion(active_set, { whichMax(gradient.cwiseAbs()) });
}

template<typename MatrixType>
bool
StrongScreening::checkKktViolationsImpl(Eigen::VectorXd& gradient,
                                        const Eigen::VectorXd& beta,
                                        const Eigen::ArrayXd& lambda_curr,
                                        std::vector<int>& working_set,
                                        const MatrixType& x,
                                        const Eigen::MatrixXd& residual,
                                        const Eigen::VectorXd& x_centers,
                                        const Eigen::VectorXd& x_scales,
                                        JitNormalization jit_normalization)
{
  // First check for violations in the strong set
  updateGradient(gradient,
                 x,
                 residual,
                 strong_set,
                 x_centers,
                 x_scales,
                 Eigen::VectorXd::Ones(x.rows()),
                 jit_normalization);

  auto violations =
    setDiff(kktCheck(gradient, beta, lambda_curr, strong_set), working_set);

  if (violations.empty()) {
    // Now check for violations in the full set
    updateGradient(gradient,
                   x,
                   residual,
                   x_centers,
                   x_scales,
                   Eigen::VectorXd::Ones(x.rows()),
                   jit_normalization);

    violations = setDiff(kktCheck(gradient, beta, lambda_curr), working_set);

    if (violations.empty()) {
      return true; // No violations found
    }
  }

  // If we found violations, update the working set
  working_set = setUnion(working_set, violations);
  return false; // Violations found
}

bool
StrongScreening::checkKktViolations(Eigen::VectorXd& gradient,
                                    const Eigen::VectorXd& beta,
                                    const Eigen::ArrayXd& lambda_curr,
                                    std::vector<int>& working_set,
                                    const Eigen::MatrixXd& x,
                                    const Eigen::MatrixXd& residual,
                                    const Eigen::VectorXd& x_centers,
                                    const Eigen::VectorXd& x_scales,
                                    JitNormalization jit_normalization)
{
  return checkKktViolationsImpl(gradient,
                                beta,
                                lambda_curr,
                                working_set,
                                x,
                                residual,
                                x_centers,
                                x_scales,
                                jit_normalization);
}

bool
StrongScreening::checkKktViolations(Eigen::VectorXd& gradient,
                                    const Eigen::VectorXd& beta,
                                    const Eigen::ArrayXd& lambda_curr,
                                    std::vector<int>& working_set,
                                    const Eigen::SparseMatrix<double>& x,
                                    const Eigen::MatrixXd& residual,
                                    const Eigen::VectorXd& x_centers,
                                    const Eigen::VectorXd& x_scales,
                                    JitNormalization jit_normalization)
{
  return checkKktViolationsImpl(gradient,
                                beta,
                                lambda_curr,
                                working_set,
                                x,
                                residual,
                                x_centers,
                                x_scales,
                                jit_normalization);
}

bool
StrongScreening::checkKktViolations(Eigen::VectorXd& gradient,
                                    const Eigen::VectorXd& beta,
                                    const Eigen::ArrayXd& lambda_curr,
                                    std::vector<int>& working_set,
                                    const Eigen::Map<Eigen::MatrixXd>& x,
                                    const Eigen::MatrixXd& residual,
                                    const Eigen::VectorXd& x_centers,
                                    const Eigen::VectorXd& x_scales,
                                    JitNormalization jit_normalization)
{
  return checkKktViolationsImpl(gradient,
                                beta,
                                lambda_curr,
                                working_set,
                                x,
                                residual,
                                x_centers,
                                x_scales,
                                jit_normalization);
}

bool
StrongScreening::checkKktViolations(
  Eigen::VectorXd& gradient,
  const Eigen::VectorXd& beta,
  const Eigen::ArrayXd& lambda_curr,
  std::vector<int>& working_set,
  const Eigen::Map<Eigen::SparseMatrix<double>>& x,
  const Eigen::MatrixXd& residual,
  const Eigen::VectorXd& x_centers,
  const Eigen::VectorXd& x_scales,
  JitNormalization jit_normalization)
{
  return checkKktViolationsImpl(gradient,
                                beta,
                                lambda_curr,
                                working_set,
                                x,
                                residual,
                                x_centers,
                                x_scales,
                                jit_normalization);
}

std::string
StrongScreening::toString() const
{
  return "strong";
}

// Factory function to create appropriate screening rule
std::unique_ptr<ScreeningRule>
createScreeningRule(const std::string& screening_type)
{
  if (screening_type == "none") {
    return std::make_unique<NoScreening>();
  } else if (screening_type == "strong") {
    return std::make_unique<StrongScreening>();
  } else {
    throw std::invalid_argument("Unknown screening type: " + screening_type);
  }
}

} // namespace slope
