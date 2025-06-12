/**
 * @file
 * @brief Implementation of screening rules for SLOPE regression optimization
 */

#include "screening.h"
#include "kkt_check.h"
#include "math.h"
#include "utils.h"
#include <Eigen/Core>
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
NoScreening::initialize(const std::vector<int>& full_set, int)
{
  return full_set;
}

std::vector<int>
NoScreening::screen(Eigen::VectorXd&,
                    const Eigen::ArrayXd&,
                    const Eigen::ArrayXd&,
                    const Eigen::VectorXd&,
                    const std::vector<int>& full_set)
{
  // No screening - use all variables
  return full_set;
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
                                JitNormalization,
                                const std::vector<int>&)
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
                                JitNormalization,
                                const std::vector<int>&)
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
                                JitNormalization,
                                const std::vector<int>&)
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
                                JitNormalization,
                                const std::vector<int>&)
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
StrongScreening::initialize(const std::vector<int>&, int alpha_max_ind)
{
  return { alpha_max_ind };
}

std::vector<int>
StrongScreening::screen(Eigen::VectorXd& gradient,
                        const Eigen::ArrayXd& lambda_curr,
                        const Eigen::ArrayXd& lambda_prev,
                        const Eigen::VectorXd& beta,
                        const std::vector<int>& full_set)
{

  if (lambda_curr(0) == 0.0) {
    // If lambda is zero, we cannot screen any features
    return full_set;
  }

  std::vector<int> active_set = activeSet(beta);
  strong_set = strongSet(gradient, lambda_curr, lambda_prev);
  strong_set = setUnion(strong_set, active_set);

  // Return working set based on active set and maximum gradient
  return setUnion(active_set, { whichMax(gradient.cwiseAbs()) });
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
                                        JitNormalization jit_normalization,
                                        const std::vector<int>& full_set)
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
                   full_set,
                   x_centers,
                   x_scales,
                   Eigen::VectorXd::Ones(x.rows()),
                   jit_normalization);

    violations =
      setDiff(kktCheck(gradient, beta, lambda_curr, full_set), working_set);

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
                                    JitNormalization jit_normalization,
                                    const std::vector<int>& full_set)
{
  return checkKktViolationsImpl(gradient,
                                beta,
                                lambda_curr,
                                working_set,
                                x,
                                residual,
                                x_centers,
                                x_scales,
                                jit_normalization,
                                full_set);
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
                                    JitNormalization jit_normalization,
                                    const std::vector<int>& full_set)
{
  return checkKktViolationsImpl(gradient,
                                beta,
                                lambda_curr,
                                working_set,
                                x,
                                residual,
                                x_centers,
                                x_scales,
                                jit_normalization,
                                full_set);
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
                                    JitNormalization jit_normalization,
                                    const std::vector<int>& full_set)
{
  return checkKktViolationsImpl(gradient,
                                beta,
                                lambda_curr,
                                working_set,
                                x,
                                residual,
                                x_centers,
                                x_scales,
                                jit_normalization,
                                full_set);
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
  JitNormalization jit_normalization,
  const std::vector<int>& full_set)
{
  return checkKktViolationsImpl(gradient,
                                beta,
                                lambda_curr,
                                working_set,
                                x,
                                residual,
                                x_centers,
                                x_scales,
                                jit_normalization,
                                full_set);
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
