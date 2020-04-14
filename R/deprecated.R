# deprecated functions ----------------------------------------------------

#' Lambda sequences for SLOPE (deprecated)
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

#' Sorted L1 solver (deprecated)
#'
#' Solves the sorted L1 penalized regression problem: given a matrix \eqn{A},
#' a vector \eqn{b}, and a decreasing vector \eqn{\lambda}, find the vector
#' \eqn{x} minimizing
#' \deqn{\frac{1}{2}\Vert Ax - b \Vert_2^2 +
#'       \sum_{i=1}^p \lambda_i |x|_{(i)}.}
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

  # Get problem dimension.
  n = ncol(A)

  # Get initial lower bound on the Lipschitz constant.
  x = with_seed(0, stats::rnorm(n))
  x = x / sqrt(sum(x^2))
  x = t(A) %*% (A %*% x)
  L = sqrt(sum(x^2))

  # Initialize parameters and iterates.
  x.init = if (is.null(initial)) rep(0,n) else initial
  t      = 1
  eta    = 2
  x      = x.init
  y      = x
  Ax     = A %*% x
  f.prev = Inf
  iter   = 0
  optimal = FALSE

  # Main loop.
  repeat {
    # Compute the gradient at f(y).
    if ((iter %% grad_iter) == 0) # Includes first iterations
      r = (A %*% y) - b
    else
      r = (Ax + ((t.prev - 1) / t) * (Ax - Ax.prev)) - b
    g = t(A) %*% r
    f = as.double(crossprod(r)) / 2

    # Increment iteration count.
    iter = iter + 1

    # Check optimality conditions.
    if ((iter %% opt_iter) == 0) {
      # Compute 'dual', check infeasibility and gap.
      gs     = sort(abs(g), decreasing=TRUE)
      ys     = sort(abs(y), decreasing=TRUE)
      infeas = max(max(cumsum(gs-lambda)),0)

      # Compute primal and dual objective.
      obj_primal =  f + as.double(crossprod(lambda,ys))
      obj_dual   = -f - as.double(crossprod(r,b))

      # Check primal-dual gap.
      if ((abs(obj_primal - obj_dual)/max(1,obj_primal) < tol_rel_gap) &&
          (infeas < tol_infeas * lambda[[1]]))
        optimal = TRUE;
    }

    # Stopping criteria.
    if (optimal || (iter >= max_iter))
      break;

    # Store copies of previous values.
    Ax.prev = Ax
    x.prev = x; f.prev = f; t.prev = t

    # Lipschitz search.
    repeat {
      # Compute prox mapping.
      x = prox(y - (1/L)*g, lambda/L)
      d = x - y

      Ax = A %*% x
      r  = Ax-b
      f  = as.double(crossprod(r))/2
      q  = f.prev + as.double(crossprod(d,g)) + (L/2)*as.double(crossprod(d))

      if (q < f*(1-1e-12))
        L = L * eta
      else
        break
    }

    # Update.
    t <- (1 + sqrt(1 + 4*t^2)) / 2
    y <- x + ((t.prev - 1) / t) * (x - x.prev)
  }
  if (!optimal)
    warning('SLOPE solver reached iteration limit')

  # Package up the results.
  structure(list(x = as.vector(y),
                 optimal = optimal,
                 iter = iter,
                 infeas = infeas,
                 obj_primal = obj_primal,
                 obj_dual = obj_dual,
                 lipschitz = L),
            class = 'SLOPE_solver.result')
}

# Evaluate an expression with the given random seed, then restore the old seed.
with_seed <- function(seed, expr) {
  seed.old <- if (exists('.Random.seed')) .Random.seed else NULL
  set.seed(seed)
  on.exit({
    if (is.null(seed.old)) {
      if (exists('.Random.seed'))
        rm(.Random.seed, envir=.GlobalEnv)
    } else {
      .Random.seed <<- seed.old
    }
  })
  expr
}

#' Prox for sorted L1 norm (deprecated)
#'
#' Compute the prox for the sorted L1 norm. That is, given a vector \eqn{x}
#' and a decreasing vector \eqn{\lambda}, compute the unique value of \eqn{y}
#' minimizing
#' \deqn{\frac{1}{2} \Vert x - y \Vert_2^2 +
#'       \sum_{i=1}^n \lambda_i |y|_{(i)}.}
#'
#' @param x input vector
#' @param lambda vector of \eqn{\lambda}'s, sorted in decreasing order
#' @param method deprecated
#'
#' @details At present, two methods for computing the sorted L1 prox are
#' supported. By default, we use a fast custom C implementation. Since SLOPE
#' can be viewed as an isotonic regression problem, the prox can also be
#' computed using the \code{isotone} package. This option is provided
#' primarily for testing.
#'
#' @export
prox_sorted_L1 <- function(x, lambda, method) {
  .Deprecated()

  sorted_l1_prox(as.matrix(x), lambda)
}
