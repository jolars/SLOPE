#' rescale
#'
#' @param x A number
#' @param y A number

rescale <- function(y,x){
  z <- data.frame(y=y,x=x)
  z <- na.omit(z)

  return(lm(x~y,data=z)$coef)
}

#' rescale_all
#'
#' @param results A number
#' @param Xmis A number

rescale_all<-function(results,Xmis){
  k = ncol(results$X)
  scales =  sapply(1:k,function(l) rescale(results$X[,l],Xmis[,l]))
  results$X = t(t(results$X)*scales[2,]+scales[1,])
  results$Sigma = results$Sigma*(scales[2,]%*%t(scales[2,]))
  results$mu = results$mu*scales[2]+scales[1,]
  results$beta = results$beta/scales[2,]
  results$intercept = -sum(scales[1,]*results$beta)
  return(results)
}

#' Adaptive Bayesian SLOPE
#'
#' @description Fit a gaussian model regularized with Adaptive Bayesian SLOPE
#' and handle missing values by Stochastic Approximation of Expected
#' Maximization (SAEM)
#'
#' @param start the initial vector of regression coefficients for the first
#' iteration. Default to the LASSO estimator obtained after
#'
#' @param Xmis words, words
#' @param Xinit words, words
#' @param Y numeric. Response variable.
#' @param a_prior,b_prior non-negative parameters of the prior Beta distribution
#' on theta.
#' @param Covmat numeric covariance matrix. Default to identity matrix.
#' @param sigma the variance of the noise. Default to 1.
#' @param FDR False Discovery Rate. Default to 0.05.
#' @param tol optimization tolerance.
#' @param known_sigma logical. Indicates whether the sigma is known
#' @param max_iter the maximal number of iterations of the optimization
#' algorithm. Default to 100.
#' @param verbose verbosity. Default to FALSE
#' @param BH logical. Indicates whether the Benjamini-Hochberg correction for
#' multiple testing should be used.
#' @param known_cov logical. Indicates whether the covariance matrix is known.
#'
#' @details \code{ABSLOPE} is the combination of SLOPE and Spike-and-Slab
#' LASSO (SSL). This approach relies on iterations of the weighted SLOPE
#' algorithm and finds the solution by minimizing
#' \deqn{
#'  G(Y,X,\beta) + \sum_{j = 1}^{p} w_j \lambda_{r(\beta,j)}|\beta_j|,
#' }
#' where \eqn{r(\beta,j) = {1,2,...,p}} is the rank of \eqn{\beta_j}  among
#' elements in \eqn{\beta} in a descending order. The weight \eqn{w_j} depends
#' on the posterior probability that a variable \eqn{X_j} is a true predictor
#' and is calculated based on the prior knowledge and on the  estimator of
#' \eqn{\beta_j}, the signal sparsity and its average strength from the previous
#' iterations.
#'
#' @references
#'
#' Jiang, W., Bogdan, M., Josse, J., Majewski, S., Miasojedow, B.,
#' Ročková, V., & TraumaBase® Group. (2021). Adaptive Bayesian SLOPE: Model
#' Selection with Incomplete Data. Journal of Computational and Graphical
#' Statistics, 1-25. \doi{10.1080/10618600.2021.1963263}
#'
#' Bogdan, M., van den Berg, E., Sabatti, C., Su, W., & Candès, E. J. (2015).
#' SLOPE -- adaptive variable selection via convex optimization. The Annals of
#' Applied Statistics, 9(3), 1103–1140. \doi{10/gfgwzt}
#'
#'Ročková, V., & George, E. I. (2018). The spike-and-slab lasso. Journal of
#'the American Statistical Association, 113(521), 431-444.
#'\doi{10.1080/01621459.2016.1260469}
#'
#' @examples
#' library(graphics)
#' library(mice)
#' X <- airquality[, c("Ozone", "Solar.R", "Wind")]
#' X <- as.matrix(X)
#' Y <- airquality$Temp
#' imp <- mice(X, m = 1, printFlag = FALSE)
#' Xinit <- as.matrix(complete(imp))
#' Covmat <- as.matrix(cov(Xinit))
#' library(glmnet)
#' obj3<-cv.glmnet(Xinit, Y,standardize=FALSE, intercept=FALSE)
#' betal<-coefficients(obj3, s='lambda.min');
#' coefs <- betal[2:length(betal)]
#' fit <- ABSLOPE(X, Y, start = coefs, Xinit = Xinit, a_prior = 0.01,
#' b_prior = 0.01, Covmat = diag(rep(1,length(start))), sigma = 1,
#' FDR = 0.05, tol = 1e-04, max_iter = 100L,
#' verbose = FALSE, BH = TRUE)
#' @export ABSLOPE
#'

ABSLOPE <- function(
  Xmis,
  Y,
  start = NULL,
  Xinit = NULL,
  a_prior,
  b_prior,
  Covmat = diag(rep(1,length(start))),
  sigma = 1,
  FDR = 0.05,
  tol = 1e-04,
  max_iter = 100L,
  verbose = FALSE,
  BH = TRUE)
{
  #TODO: checkmate


  # if Covmat is null -> known_cov = FALSE
  known_cov <- !is.null(Covmat)
  # if sigma is null -> known_sigma = FALSE
  known_sigma <- !is.null(sigma)

  #TODO: if Xinit is null -> mice imputation
  if (is.null(Xinit)) {
    #imp = mice(Xtrauma,m=1, printFlag = FALSE)
    #Xfull.sim  = complete(imp)
    #save(Xfull.sim, file="imputed_mice.Rdata")

  }

  #TODO: if start is null -> LASSO gives starting coefficients
  if (is.null(start)) {
    invisible()
  }

  out <- SLOBE_ADMM_approx_missing(start, Xmis, Xinit, Y, a_prior, b_prior,
                                   Covmat, sigma, FDR , tol, known_sigma,
                                   max_iter, verbose , BH , known_cov )
  return(list(out, Xmis))
  #TODO: fix this!
  #return(rescale_all(out,Xmis))
}

# TODO:
#
# 1. If start == NULL then calculate LASSO based on simple imputation with averaged values
#
# 2. If Xinit == NULL then use default imputation with means
# Do we want to use the parameter Xinit or the initial imputation method or both functionalities?
#
# 3. What default values for a_prior and b_prior should we choose?
#
# 4. If sigma is NULL then assume that known_sigma is FALSE. Known_sigma is redundant
# (the same as in the case of covmat)
#
# 5. verbose is in polish (iteracja)
#
# 6. vector lambda is hardcoded
#
# 7. examples




