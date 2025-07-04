% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpolate_coefficients.R
\name{interpolateCoefficients}
\alias{interpolateCoefficients}
\title{Interpolate coefficients}
\usage{
interpolateCoefficients(beta, intercepts, interpolation_list)
}
\arguments{
\item{beta}{coefficients}

\item{intercepts}{intercepts}

\item{interpolation_list}{a list generated from \code{\link[=interpolatePenalty]{interpolatePenalty()}}}
}
\value{
A matrix (or list of matrices) with new coefficients based
on linearly interpolating from new and old lambda values.
}
\description{
Interpolate coefficients
}
\keyword{internal}
