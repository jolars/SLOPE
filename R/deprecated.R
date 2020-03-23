#' Lambda sequences for SLOPE
#'
#' Computes \eqn{\lambda} sequences for SLOPE according to several pre-defined methods.
#'
#' @param n number of observations
#' @param p number of variables
#' @param fdr target False Discovery Rate (FDR)
#' @param method method to use for computing \eqn{\lambda} (see Details)
#'
#' @details The following methods for computing \eqn{\lambda} are supported:
#' \itemize{
#'  \item \code{bhq}: Computes sequence inspired by Benjamini-Hochberg (BHq)
#'    procedure
#'  \item \code{gaussian}: Computes modified BHq sequence inspired by
#'    Gaussian designs
#' }
#'
#' @rdname lambda
#' @export
create_lambda <- function(n, p, fdr = 0.20, method = c('bhq','gaussian')) {

  .Deprecated()

  impl = switch(match.arg(method),
                bhq = create_lambda_bhq,
                gaussian = create_lambda_gaussian_truncated)
  impl(n, p, fdr)
}

create_lambda_bhq <- function(n, p, fdr) {
  q = (1:p) * fdr / (2*p)
  stats::qnorm(1 - q)
}

create_lambda_gaussian <- function(n, p, fdr) {
  w <- function(k) 1 / max(1, n - k - 1)
  lambda.bhq = create_lambda_bhq(n, p, fdr)
  lambda = rep(0,p)
  lambda[1] <- lambda.bhq[1]
  if (p >= 2) {
    sum_sq <- 0
    for (i in 2:p) {
      sum_sq <- sum_sq + lambda[i-1]^2
      lambda[i] <- lambda.bhq[i] * sqrt(1 + w(i-1) * sum_sq)
    }
  }
  lambda
}

create_lambda_gaussian_truncated <- function(n, p, fdr) {
  lambda = create_lambda_gaussian(n, p, fdr)
  k = which.min(lambda)
  lambda[k:p] <- lambda[k]
  lambda
}


#' Prox for sorted L1 norm
#'
#' Compute the prox for the sorted L1 norm. That is, given a vector \eqn{x}
#' and a decreasing vector \eqn{\lambda}, compute the unique value of \eqn{y}
#' minimizing
#' \deqn{\frac{1}{2} \Vert x - y \Vert_2^2 +
#'       \sum_{i=1}^n \lambda_i |x|_{(i)}.}
#'
#' @param x input vector
#' @param lambda vector of \eqn{\lambda}'s, sorted in decreasing order
#' @param method underlying prox implementation, either 'c' or 'isotone'
#'  (see Details)
#'
#' @details At present, two methods for computing the sorted L1 prox are
#' supported. By default, we use a fast custom C implementation. Since SLOPE
#' can be viewed as an isotonic regression problem, the prox can also be
#' computed using the \code{isotone} package. This option is provided
#' primarily for testing.
#'
#' @export
prox_sorted_L1 <- function(x, lambda, method = c('c','isotone')) {
  .Deprecated()

  sorted_l1_prox(x, lambda);
}

#' Sorted L1 solver (deprecated)
#'
#' This function is deprecated. Please call [SLOPE()] directly instead.
#'
#' @param A an \eqn{n}-by-\eqn{p} matrix
#' @param b vector of length \eqn{n}
#' @param lambda vector of length \eqn{p}, sorted in decreasing order
#' @param initial initial guess for \eqn{x}
#' @param prox function that computes the sorted L1 prox
#' @param max_iter maximum number of iterations in the gradient descent
#' @param grad_iter number of iterations between gradient updates
#' @param opt_iter number of iterations between checks for optimality
#' @param tol_infeas tolerance for infeasibility
#' @param tol_rel_gap tolerance for relative gap between primal and dual
#'  problems
#'
#' @return An object of class \code{SLOPE_solver.result}. This object is a list
#'  containing at least the following components:
#'  \item{x}{solution vector \eqn{x}}
#'  \item{optimal}{logical: whether the solution is optimal}
#'  \item{iter}{number of iterations}
#'
#' @details This optimization problem is convex and is solved using an
#' accelerated proximal gradient descent method.
#'
#' @export
# Adapted from Adlas.R in original SLOPE distribution.
SLOPE_solver <- function(A, b, lambda, initial = NULL, prox = prox_sorted_L1,
                         max_iter = 10000, grad_iter = 20, opt_iter = 1,
                         tol_infeas = 1e-6, tol_rel_gap = 1e-6) {
  .Deprecated("SLOPE")

  SLOPE(A,
        b,
        lambda = lambda,
        max_passes = max_iter,
        tol_infeas = tol_infeas,
        tol_rel_gap = tol_rel_gap)
}
