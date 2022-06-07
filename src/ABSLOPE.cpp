/*
 * Copyright 2019, Szymon Majewski, Błażej Miasojedow
 *
 * This file is part of SLOBE-Rcpp Toolbox
 *
 *   The SLOBE-Rcpp Toolbox is free software: you can redistribute it
 *   and/or  modify it under the terms of the GNU General Public License
 *   as published by the Free Software Foundation, either version 3 of
 *   the License, or (at your option) any later version.
 *
 *   The SLOPE Toolbox is distributed in the hope that it will
 *   be useful, but WITHOUT ANY WARRANTY; without even the implied
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *   See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the SLOPE Toolbox. If not, see
 *   <http://www.gnu.org/licenses/>.
 */

#include <RcppArmadillo.h>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <stdlib.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

template<typename T>

/*
 * Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes
 *
 * The following function is copied from proxSortedL1.c file
 * which is part of SLOPE Toolbox.
 *
 *   The SLOPE Toolbox is free software: you can redistribute it
 *   and/or  modify it under the terms of the GNU General Public License
 *   as published by the Free Software Foundation, either version 3 of
 *   the License, or (at your option) any later version.
 *
 *   The SLOPE Toolbox is distributed in the hope that it will
 *   be useful, but WITHOUT ANY WARRANTY; without even the implied
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *   See the GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the SLOPE Toolbox. If not, see
 *   <http://www.gnu.org/licenses/>.
 */

int
evaluateProx(double* y, double* lambda, double* x, size_t n, int* order);

/* ----------------------------------------------------------------------- */
int
evaluateProx(double* y, double* lambda, double* x, size_t n, int* order)
{
  double d;
  double* s = NULL;
  double* w = NULL;
  size_t* idx_i = NULL;
  size_t* idx_j = NULL;
  size_t i, j, k;
  int result = 0;

  /* Allocate memory */
  s = (double*)malloc(sizeof(double) * n);
  w = (double*)malloc(sizeof(double) * n);
  idx_i = (size_t*)malloc(sizeof(size_t) * n);
  idx_j = (size_t*)malloc(sizeof(size_t) * n);

  if ((s != NULL) && (w != NULL) && (idx_i != NULL) && (idx_j != NULL)) {
    k = 0;
    for (i = 0; i < n; i++) {
      idx_i[k] = i;
      idx_j[k] = i;
      s[k] = y[i] - lambda[i];
      w[k] = s[k];

      while ((k > 0) && (w[k - 1] <= w[k])) {
        k--;
        idx_j[k] = i;
        s[k] += s[k + 1];
        w[k] = s[k] / (i - idx_i[k] + 1);
      }

      k++;
    }

    if (order == NULL) {
      for (j = 0; j < k; j++) {
        d = w[j];
        if (d < 0)
          d = 0;

        for (i = idx_i[j]; i <= idx_j[j]; i++) {
          x[i] = d;
        }
      }
    } else {
      for (j = 0; j < k; j++) {
        d = w[j];
        if (d < 0)
          d = 0;
        for (i = idx_i[j]; i <= idx_j[j]; i++) {
          x[order[i]] = d;
        }
      }
    }
  } else {
    result = -1;
  }

  /* Deallocate memory */
  if (s != NULL)
    free(s);
  if (w != NULL)
    free(w);
  if (idx_i != NULL)
    free(idx_i);
  if (idx_j != NULL)
    free(idx_j);

  return result;
}

// Comparator class for argsort function
class CompareByNumericVectorValues
{
private:
  const NumericVector* _values;

public:
  CompareByNumericVectorValues(const NumericVector* values)
  {
    _values = values;
  }
  bool operator()(const int& a, const int& b)
  {
    return ((*_values)[a] > (*_values)[b]);
  }
};

// Writes down in IntegerVector ord sequnce of indexes,
// such that w[ord[i]] >= w[ord[j]] whenever i <= j
// for the given NumericVector w.
void
argsort(const NumericVector& w, IntegerVector ord)
{
  std::iota(ord.begin(), ord.end(), 0);
  CompareByNumericVectorValues comp = CompareByNumericVectorValues(&w);
  std::sort(ord.begin(), ord.end(), comp);
}

// Computes proximal step of SLOPE(lambda) norm from point y
// NumericVector y is not assumed to be sorted, sorting is performed within
// function
NumericVector
prox_sorted_L1_C(NumericVector y, NumericVector lambda)
{
  int n = y.size();
  NumericVector x(n);
  IntegerVector order(n);
  argsort(abs(y), order);
  IntegerVector sign_y = sign(y);
  y = abs(y);
  y.sort(true);
  evaluateProx(y.begin(), lambda.begin(), x.begin(), n, NULL);
  NumericVector res(n);
  for (int k = 0; k < n; k++) {
    res[order[k]] = sign_y[order[k]] * x[k];
  }
  return res;
}

// Creates a vector of weights for SLOPE for a given p and FDR
// Writes down the vector in the passed NumericVector lam
void
createLambda(NumericVector& lam, int p, double fdr, bool bh)
{
  NumericVector h(p);

  if (bh) {
    for (double i = 0.0; i < h.size(); ++i) {
      h[i] = 1 - (fdr * (i + 1) / (2 * p));
    }
  } else {
    for (double i = 0.0; i < h.size(); ++i) {
      h[i] = 1 - (fdr * (1) / (2 * p));
    }
  }

  lam = qnorm(h);
}

// Expectation of truncated gamma distribution
double
exTruncGamma(double a, double b)
{
  double c = exp(Rf_pgamma(1.0, a + 1, 1.0 / b, 1, 1) -
                 Rf_pgamma(1.0, a, 1.0 / b, 1, 1));
  c /= b;
  c *= a;

  return c;
}

// Compute SLOBE estimator for a linear model using ADMM
arma::vec
slopeADMM(const arma::mat& x,
          const vec& y,
          NumericVector& lambda,
          const int& p,
          const double& rho,
          int max_iter = 500,
          double tol_inf = 1e-08)
{
  // Precompute M = (X^TX + rho I)^{-1}
  // and MXtY = M * X^T * Y for proximal steps of quadratic part
  arma::mat xtx = x.t() * x;
  for (int i = 0; i < p; ++i) {
    xtx(i, i) += rho;
  }
  xtx = xtx.i();
  arma::vec xtx_xty = xtx * (x.t() * y);
  NumericVector lam_seq_rho = lambda / rho;

  // Prepare variables before starting ADMM loop
  int i = 0;
  arma::vec beta = zeros(p);
  arma::vec z = zeros(p);
  arma::vec u = zeros(p);
  NumericVector z_new = NumericVector(p);
  arma::vec z_new_arma = zeros(p);
  NumericVector beta_plus_u(p);
  double dual_feas, primal_feas;

  // ADMM loop
  while (i < max_iter) {
    beta = xtx_xty + xtx * (rho * (z - u));
    beta_plus_u = as<NumericVector>(wrap(beta + u));
    z_new = prox_sorted_L1_C(beta_plus_u, lam_seq_rho);
    z_new_arma = as<arma::vec>(z_new);
    u += (beta - z_new_arma);

    dual_feas = arma::norm(rho * (z_new_arma - z));
    primal_feas = arma::norm(z_new_arma - beta);

    z = z_new_arma;

    if (primal_feas < tol_inf && dual_feas < tol_inf) {
      i = max_iter;
    }

    ++i;
  }

  return z;
}

// Scaling matrix X by weight vector w
void
divXbyW(mat& x_div_w,
        const mat& X,
        const vec& w_vec,
        const int& n,
        const int& p)
{

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      x_div_w(i, j) = X(i, j) / w_vec(j);
    }
  }
}

// Missing data imputation procedure that imputes means of non-missing values
// in a column
void
imputeMean(mat& x, const int& n, const int& p)
{

  double colmean;
  int non_na_rows;
  for (int c_num = 0; c_num < p; c_num++) {
    colmean = 0.0;
    non_na_rows = 0;

    for (int r_num = 0; r_num < n; r_num++) {
      if (arma::is_finite(x(r_num, c_num))) {
        colmean += x(r_num, c_num);
        non_na_rows += 1;
      }
    }
    colmean /= non_na_rows;
    for (int r_num = 0; r_num < n; r_num++) {
      if (!arma::is_finite(x(r_num, c_num))) {
        x(r_num, c_num) = colmean;
      }
    }
  }
}

void
linShrinkCov(mat& x, arma::mat& s_mat, const int& n, const int& p)
{
  rowvec means = sum(x) / n;
  s_mat = x.t() * x;
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      s_mat(i, j) -= n * means[i] * means[j];
      s_mat(i, j) /= n;
    }
  }
  double m = trace(s_mat) / p;
  double d2 = norm(s_mat, "fro");
  double b_bar2 = pow(d2, 2);
  d2 = (b_bar2 - p * pow(m, 2)) / p;
  b_bar2 *= n;
  double prod;
  double sum_prod;
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      sum_prod = 0.0;
      for (int k = 0; k < n; ++k) {
        prod = (x(k, i) - means[i]) * (x(k, j) - means[j]);
        b_bar2 += prod * prod;
        sum_prod += prod;
      }
      b_bar2 -= 2 * s_mat(i, j) * sum_prod;
    }
  }
  b_bar2 /= (p * pow(n - 1, 2));
  double b2 = b_bar2 < d2 ? b_bar2 : d2;
  double a2 = d2 - b2;
  s_mat *= (a2 / d2);
  m *= (b2 / d2);
  for (int i = 0; i < p; ++i) {
    s_mat(i, i) += m;
  }
}

void
imputeRowAdvance(const vec& beta,
                 mat& X,
                 const vec& Y,
                 const mat& S,
                 const double& sigma_sq,
                 const LogicalMatrix& XisFin,
                 const int& n,
                 const int& p,
                 const int& row,
                 const std::vector<int> nanCols,
                 const vec& m,
                 const vec& tau_sq)
{

  int l = nanCols.size();
  int s, t;
  arma::mat A = zeros(l, l);
  arma::vec u = zeros(l);
  double r = Y[row];
  int u_ind = 0;
  for (int i = 0; i < p; ++i) {
    if (!XisFin(row, i)) {
      for (int j = 0; j < p; ++j) {
        if (XisFin(row, j)) {
          u[u_ind] += X(row, j) * S(j, i);
        }
      }
      u_ind += 1;
    } else {
      r -= beta[i] * X(row, i);
    }
  }

  for (int i = 0; i < l; ++i) {
    for (int j = 0; j < l; ++j) {
      if (i == j) {
        A(i, j) = 1.0;
      } else {
        s = nanCols[i];
        t = nanCols[j];
        A(i, j) = (beta[s] * beta[t] / sigma_sq + S(s, t)) / tau_sq[s];
      }
    }
  }

  arma::vec b = zeros(l);
  for (int i = 0; i < l; ++i) {
    t = nanCols[i];
    b[i] = ((r * beta[t]) / sigma_sq + m[t] - u[i]) / tau_sq[t];
  }

  arma::vec sol = solve(A, b);
  for (int i = 0; i < l; ++i) {
    t = nanCols[i];
    X(row, t) = sol[i];
  }
}

// Imputation procedure for SLOBE with missing values
void
imputeAdvance(const arma::vec& beta,
              mat& x,
              const arma::vec& y,
              const mat& s_mat,
              const double& sigma_sq,
              const int& n,
              const int& p,
              const vec& mu,
              const LogicalMatrix& x_is_fin,
              const std::vector<int> any_nan_x_rows,
              const std::vector<std::vector<int>> nan_ind_in_rows)
{

  arma::vec tau_sq = zeros(p);
  tau_sq = square(beta) / sigma_sq;
  for (int i = 0; i < p; ++i) {
    tau_sq[i] = tau_sq[i] + s_mat(i, i);
  }

  arma::vec m = zeros(p);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      m[i] += mu[j] * s_mat(i, j);
    }
  }

  for (size_t i = 0; i < any_nan_x_rows.size(); ++i) {
    imputeRowAdvance(beta,
                     x,
                     y,
                     s_mat,
                     sigma_sq,
                     x_is_fin,
                     n,
                     p,
                     any_nan_x_rows[i],
                     nan_ind_in_rows[i],
                     m,
                     tau_sq);
  }
}

// A subroutine that updates gamma by using the mean lambda when several beta_i
// have the same magnitude
void
updateGammaMean(const NumericVector& abs_beta_ord,
                const NumericVector& lambda,
                const int& p,
                NumericVector& gamma_h)
{

  double bl_sum = lambda[0] * abs_beta_ord[0];
  double bl_mean;
  int equals_in_sequence = 1;
  for (int i = 1; i < p; ++i) {
    if (abs_beta_ord[i] == abs_beta_ord[i - 1]) {
      bl_sum += lambda[i] * abs_beta_ord[i];
      equals_in_sequence += 1;
    } else {
      bl_mean = bl_sum / equals_in_sequence;
      for (int j = i - equals_in_sequence; j < i; j++) {
        gamma_h[j] = bl_mean;
      }
      bl_sum = lambda[i] * abs_beta_ord[i];
      equals_in_sequence = 1;
    }
  }
  bl_mean = bl_sum / equals_in_sequence;
  for (int j = p - equals_in_sequence; j < p; j++) {
    gamma_h[j] = bl_mean;
  }
}

// [[Rcpp::export]]
void
centerAndScale(arma::mat& x, const int& n, const int& p)
{
  for (int i = 0; i < p; ++i) {
    x.col(i) -= sum(x.col(i)) / n;
  }
  x = normalise(x);
}

// [[Rcpp::export]]
List
SLOBE_ADMM_approx_missing(NumericVector beta_start,
                          arma::mat x_miss,
                          NumericMatrix x_init,
                          arma::vec y,
                          double a_prior,
                          double b_prior,
                          arma::mat covmat,
                          double sigma = 1.0,
                          double fdr = 0.05,
                          double tol = 1e-04,
                          bool known_sigma = false,
                          int max_iter = 100,
                          bool verbose = false,
                          bool bh = true,
                          bool known_cov = false)
{

  // Initialize variables
  int p = beta_start.length();
  int n = y.size();
  NumericVector beta = clone(beta_start);
  NumericVector beta_new(p);
  arma::vec beta_arma = as<vec>(beta);
  NumericVector w(p, 1.0);
  arma::vec w_vec = ones(p);
  NumericVector wbeta(p);
  NumericVector gamma(p);
  NumericVector gamma_h(p);
  NumericVector b_sum_h(p);
  NumericVector lambda_sigma(p);
  IntegerVector order(p);
  arma::mat x_div_w = zeros(n, p);
  double error = 0.0;
  double swlambda = 0.0;
  double RSS = 0.0;
  double sigma_sq = 1.0;

  // First imputation
  bool any_nan_in_row;

  LogicalMatrix x_is_fin(n, p);
  std::vector<int> any_nan_x_rows;
  std::vector<std::vector<int>> nan_indices_in_row;
  for (int i = 0; i < n; ++i) {
    any_nan_in_row = false;
    std::vector<int> nan_ind;
    for (int j = 0; j < p; ++j) {
      x_is_fin(i, j) = arma::is_finite(x_miss(i, j));
      if (!x_is_fin(i, j)) {
        nan_ind.push_back(j);
        any_nan_in_row = true;
      }
    }
    if (any_nan_in_row) {
      any_nan_x_rows.push_back(i);
      nan_indices_in_row.push_back(nan_ind);
    }
  }

  NumericMatrix temp_x = clone(x_init);
  arma::mat x = as<mat>(temp_x);
  centerAndScale(x, n, p);
  arma::mat sigma_mat = zeros(p, p);
  if (!known_cov) {
    linShrinkCov(x, sigma_mat, n, p);
  } else {
    sigma_mat = covmat;
  }
  arma::mat S = inv_sympd(sigma_mat);
  arma::vec mu = zeros(p);
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < n; ++j) {
      mu(i) += x(j, i);
    }
    mu[i] /= n;
  }

  // Compute vector lambda based on BH procedure
  NumericVector lambda(p);
  createLambda(lambda, p, fdr, bh);

  // Initialize sigma, c, theta, gamma
  double sstart = sum(beta_start != 0);
  double pstart = 1.0;
  if (sstart > 0) {
    pstart = sstart;
  }
  if (!known_sigma) {
    sigma = sqrt(sum(pow(x * beta_arma - y, 2)) / (n - pstart));
  }
  lambda_sigma = lambda * sigma;
  argsort(beta, order);
  double c = 0.0;
  if (sstart > 0) {
    double h = (sstart + 1) / (sum(abs(beta))) * (sigma / lambda[p - 1]);
    c = (h < 1.0) ? h : 1.0;
  } else {
    c = 1.0;
  }

  double theta = (sstart + a_prior) / (a_prior + b_prior + p);

  updateGammaMean(abs(beta[order]), lambda, p, gamma_h);
  gamma_h = gamma_h * (c - 1) / sigma;

  gamma_h = (theta * c) / (theta * c + (1 - theta) * exp(gamma_h));
  for (int i = 0; i < p; ++i) {
    gamma[order[i]] = gamma_h[i];
  }

  // Start main loop
  bool converged = false;
  int iter = 0;
  while (iter < max_iter) {
    if (verbose) {
      Rcout << "Iteration: " << iter << "/" << max_iter << "\n";
    }
    w = 1.0 - (1.0 - c) * gamma;

    // Compute rewieghted SLOPE estimator using computed weights and sigma
    w_vec = as<vec>(w);
    divXbyW(x_div_w, x, w_vec, n, p);
    beta_arma = slopeADMM(x_div_w, y, lambda_sigma, p, 1.0);

    wbeta = as<NumericVector>(wrap(beta_arma));
    wbeta = abs(wbeta);
    argsort(wbeta, order);

    for (int i = 0; i < p; ++i) {
      beta_arma[i] /= w_vec[i];
    }
    beta_new = as<NumericVector>(wrap(beta_arma));

    // For the version with unknown sigma, estimate it
    if (!known_sigma) {
      RSS = sum(pow((x * beta_arma - y), 2));
      swlambda = sum(wbeta.sort(true) * lambda);
      sigma = (swlambda + sqrt(pow(swlambda, 2.0) + 4 * (n + 2) * RSS)) /
              (2 * (n + 2));
    }
    lambda_sigma = lambda * sigma;
    sigma_sq = sigma * sigma;

    // Estimate covariance matrix
    if (!known_cov) {
      linShrinkCov(x, sigma_mat, n, p);

      S = inv_sympd(sigma_mat);
    }
    mu = zeros(p);
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < n; ++j) {
        mu(i) += x(j, i);
      }
      mu[i] /= n;
    }

    imputeAdvance(beta_new,
                  x,
                  y,
                  S,
                  sigma_sq,
                  n,
                  p,
                  mu,
                  x_is_fin,
                  any_nan_x_rows,
                  nan_indices_in_row);

    centerAndScale(x, n, p);

    // compute new gamma
    updateGammaMean(abs(beta[order]), lambda, p, gamma_h);
    gamma_h = gamma_h * (c - 1) / sigma;

    gamma_h = (theta * c) / (theta * c + (1 - theta) * exp(gamma_h));

    for (int i = 0; i < p; ++i) {
      gamma[order[i]] = gamma_h[i];
    }

    // update c
    double sum_gamma = sum(gamma);
    b_sum_h = gamma[order];
    b_sum_h = b_sum_h * abs(beta[order]);
    b_sum_h = b_sum_h * lambda;
    double b_sum = sum(b_sum_h) / sigma;

    if (sum_gamma > 0) {
      if (b_sum > 0) {
        c = exTruncGamma(sum_gamma, b_sum);
      } else {
        c = (sum_gamma + 1) / (sum_gamma + 2);
      }
    } else {
      c = 0.5;
    }

    // update theta
    theta = (sum_gamma + a_prior) / (p + a_prior + b_prior);

    // Check stop condition
    error = sum(abs(beta - beta_new));
    if (error < tol) {
      iter = max_iter;
      converged = true;
    }
    if (verbose) {
      Rcout << "Error =  " << error << " sigma = " << sigma
            << ", theta = " << theta << ", c = " << c << "\n";
    }
    std::copy(beta_new.begin(), beta_new.end(), beta.begin());
    ++iter;
  }

  return List::create(Named("coefficients") = beta,
                      Named("sigma") = sigma,
                      Named("theta") = theta,
                      Named("c") = c,
                      Named("w") = w,
                      Named("converged") = converged,
                      Named("X") = x,
                      Named("Sigma") = sigma_mat,
                      Named("mu") = mu,
                      Named("lambda") = lambda);
}
