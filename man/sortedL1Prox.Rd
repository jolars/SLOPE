% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sorted_l1_prox.R
\name{sortedL1Prox}
\alias{sortedL1Prox}
\title{Sorted L1 Proximal Operator}
\source{
M. Bogdan, E. van den Berg, Chiara Sabatti, Weijie Su, and Emmanuel J.
Candès, “SLOPE – adaptive variable selection via convex optimization,” Ann
Appl Stat, vol. 9, no. 3, pp. 1103–1140, 2015.
}
\usage{
sortedL1Prox(x, lambda, method)
}
\arguments{
\item{x}{A vector. In SLOPE, this is the vector of coefficients.}

\item{lambda}{A non-negative and decreasing sequence
of weights for the Sorted L1 Norm. Needs to be the same length as
\code{x}.}

\item{method}{DEPRECATED}
}
\value{
An evaluation of the proximal operator at \code{x} and \code{lambda}.
}
\description{
The proximal operator for the Sorted L1 Norm, which is the penalty function
in SLOPE. It solves the problem
\deqn{
  \arg\,\min_x
    \Big(J(x, \lambda) + \frac{1}{2} ||x - v||_2^2\Big)
}{
  argmin_x (J(x, \lambda) + 0.5||x - v||_2^2)
}
where \eqn{J(x, \lambda)} is the Sorted L1 Norm.
}
