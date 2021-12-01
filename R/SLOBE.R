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

#' SLOBEC
#' @param start words, words
#' @param Xmis A number
#'
#' @export SLOBEC
#'

SLOBEC <- function(start, Xmis, Xinit, Y, a_prior, b_prior, Covmat = diag(rep(1,length(start))),sigma = 1,
                   FDR = 0.05, tol = 1e-04, known_sigma = FALSE, max_iter = 100L,
                   verbose = FALSE, BH = TRUE, known_cov = FALSE){

  out <- SLOBE_ADMM_approx_missing(start, Xmis, Xinit, Y, a_prior, b_prior, Covmat, sigma,
                                   FDR , tol, known_sigma, max_iter,
                                   verbose , BH , known_cov )

  return(rescale_all(out,Xmis))
}
