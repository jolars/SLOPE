#Formated with tidyverse
library(SLOPE)
library(truncdist)
library(nlshrink)
library(MASS)
library(glmnet)
library(missMDA)
library(mice)
# library(matlib)

ABSLOPE <- function(X, y, lambda, a = 1, b = 10, beta.start = NA, maxit = 300, case = "MCAR", seed = NA, print_iter = FALSE, tol_em = 1e-6, impute = "PCA", sigma.known = NA, sigma.init = NA, scale = FALSE, method_na = "lineq") {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  # missing pattern
  rindic <- as.matrix(is.na(X))
  if (sum(rindic) > 0) { # missing data exist
    whichcolmissing <- (1:ncol(rindic))[apply(rindic, 2, sum) > 0]
    missingcols <- length(whichcolmissing)
  }
  if (sum(rindic) == 0) {
    missingcols <- 0
  } # no missingness

  p <- ncol(X)

  # delete rows completely missing
  if (missingcols != 0) {
    if (any(apply(is.na(X), 1, sum) == p)) {
      i_allNA <- which(apply(is.na(X), 1, sum) == p)
      X <- X[-i_allNA, ]
      y <- y[-i_allNA]
    }
    if (any((is.na(y)) == TRUE)) {
      i_YNA <- which(is.na(y) == TRUE)
      X <- X[-i_YNA, ]
      y <- y[-i_YNA]
    }
  }

  n <- length(y)

  if (missingcols != 0) {
    # impute by PCA
    if (impute == "PCA") {
      X.sim <- imputePCA(X)$completeObs
    } else { # impute by mean
      X.mean <- X
      for (i in 1:ncol(X.mean)) {
        X.mean[is.na(X.mean[, i]), i] <- mean(X.mean[, i], na.rm = TRUE)
      }
      X.sim <- X.mean
    }
  } else {
    X.sim <- X
  } # no missingness


  # define functions to perform scaling
  scale_mean.w <- function(V, weight) {
    res <- sum(V * weight, na.rm = TRUE) / sum(weight[!is.na(V)])
  }
  scale_std.w <- function(V, weight) {
    res <- sqrt(sum(V^2 * weight, na.rm = TRUE) / sum(weight[!is.na(V)]))
  }
  # scaling
  if (scale) {
    row.w <- rep(1, nrow(X.sim)) / nrow(X.sim)
    mean.w <- apply(X.sim, 2, scale_mean.w, row.w)
    X.sim <- t(t(X.sim) - mean.w)
    std.w <- apply(X.sim, 2, scale_std.w, row.w) * sqrt(n)
    X.sim <- t(t(X.sim) / std.w)
  }



  ## ------------------------------------
  # initialization

  # beta
  if (is.na(beta.start)) { # initialization not given
    # LASSO cv
    objstart <- cv.glmnet(X.sim, y, standardize = FALSE, intercept = FALSE)
    coeff <- coef(objstart, s = "lambda.min")
    beta <- coeff[2:(p + 1), 1]
  } else {
    beta <- beta.start
  } # initialization given

  # sigma
  if (is.na(sigma.known)) { # real value of sigma is unknown
    if (is.na(sigma.init)) { # initialization not given
      # sigma = sd(y - X.sim %*% beta)
      nnon <- length(which(abs(beta) > 0))
      sigma <- sqrt(sum((y - X.sim %*% beta)^2) / (n - nnon))
    } else {
      sigma <- sigma.init
    } # initialization given
  } else {
    sigma <- sigma.known
  } # real value of sigma is known

  # theta
  theta <- (sum(beta != 0) + a) / (a + b + p)


  # rank
  rk <- p - rank(abs(beta), ties.method = "average") + 1

  # lambda
  lambda_sigma <- lambda * sigma
  lambda_sigma_inv <- lambda / sigma

  # c
  c <- min((sum(abs(beta) > 0) + 1) / sum(abs(beta[beta != 0])) / lambda_sigma_inv[p], 1)
  if (!is.finite(c)) {
    c <- 1
  }

  # gamma
  gamma <- (abs(beta) > 0) * 1
  pi.binom <- 1 / (1 + (1 - theta) / theta / c * exp(-abs(beta) * lambda_sigma_inv[rk] * (1 - c)))
  # gamma = rbinom(1,1,pi.binom)

  # if NA -> mu & Sigma
  if (missingcols != 0) {
    mu <- apply(X.sim, 2, mean)
    Sigma <- linshrink_cov(X.sim)
  }



  ## ------------------------------------

  # Iteration of SAEM
  # store the estimation of each iteration
  seqbeta <- seqgamma <- matrix(NA, nrow = ncol(X), ncol = (maxit + 1))
  seqsigma <- seqc <- seqtheta <- matrix(NA, nrow = 1, ncol = (maxit + 1))
  seqbeta[, 1] <- beta
  seqsigma[, 1] <- sigma
  seqgamma[, 1] <- gamma
  seqc[, 1] <- c
  seqtheta[, 1] <- theta

  cstop <- 1
  t <- 0
  while ((cstop > tol_em) * (t < maxit) | (t < 20)) {
    t <- t + 1

    if ((print_iter == TRUE) & (t %% 100 == 0)) {
      cat(sprintf("iteration = %i ", t))
      cat(sprintf("beta ="), beta)
      cat(sprintf(" sigma ="), sigma, "\n")
      cat(sprintf("Distance from last iter ="), cstop, "\n")
      # cat(sprintf(' c ='),c,'\n')
      # cat(sprintf(' pi.binom ='),pi.binom)
      # cat(sprintf(' gamma ='),gamma, '\n')
    }

    # step size
    if (t <= 20) {
      eta <- 1
    } else {
      eta <- 1 / (t - 20)
    }

    # Simulation step
    # 1) gamma
    W <- gamma * c + (rep(1, p) - gamma)
    pi.binom <- 1 / (1 + (1 - theta) / theta / c * exp(-abs(beta) * lambda_sigma_inv[rk] * (1 - c)))
    gamma <- rbinom(p, 1, pi.binom)

    # 2) c
    W <- gamma * c + (rep(1, p) - gamma)
    a.gamma <- 1 + sum(gamma == 1)
    b.gamma <- sum(abs(beta) * lambda_sigma_inv[rk] * (gamma == 1))
    # value0 = pgamma(0, shape=a.gamma, rate=b.gamma)
    # value1 = pgamma(1, shape=a.gamma, rate=b.gamma)
    # if(value0 != value1){c = rtrunc(1, "gamma", 0, 1, shape=a.gamma, rate=b.gamma)
    # }else{u = runif(1, min = 0, max = 1); c = u^(1/a.gamma)}
    if (a.gamma > 1) {
      if (b.gamma > 0) {
        c <- rtrunc(1, "gamma", 0, 1, shape = a.gamma, rate = b.gamma)
      } else {
        c <- rbeta(1, shape1 = a.gamma, shape2 = 1)
      }
    }
    else {
      c <- runif(1, 0, 1)
    }

    # 3) if NA -> Xmis
    if (missingcols != 0) {
      if (method_na == "MH") {
        S.inv <- solve(Sigma)
        for (i in (1:n)) {
          yi <- y[i]
          jna <- which(is.na(X[i, ]))
          njna <- length(jna)
          if (njna > 0) {
            xi <- X.sim[i, ]
            Oi <- Sigma[jna, jna]
            mi <- mu[jna]
            if (njna < p) {
              jobs <- setdiff(1:p, jna)
              mi <- mi + Sigma[jna, jobs] %*% solve(Sigma[jobs, jobs]) %*% (xi[jobs] - mu[jobs])
              Oi <- Oi - Sigma[jna, jobs] %*% solve(Sigma[jobs, jobs]) %*% Sigma[jobs, jna]
            }
            nmcmc <- 20
            xina <- xi[jna]
            for (m in (1:nmcmc)) {
              xina.c <- mvrnorm(n = 1, mu = mi, Sigma = Oi)
              xi[jna] <- xina.c
              alpha <- dnorm(yi, mean = xi %*% beta, sd = sigma, log = FALSE) / dnorm(yi, mean = xi[jobs] %*% beta[jobs] + xina %*% beta[jna], sd = sigma, log = FALSE)
              if (runif(1) < alpha) {
                xina <- xina.c
              }
            }
            X.sim[i, jna] <- xina
          }
        }
      }

      if (method_na == "lineq") {
        S.inv <- solve(Sigma)
        m <- S.inv %*% mu
        tau <- sqrt(diag(S.inv) + (beta / sigma)^2)
        for (i in (1:n)) {
          yi <- y[i]
          jna <- which(is.na(X[i, ]))
          njna <- length(jna)
          if (njna > 0) {
            xo <- X[i, -jna]
            betai <- beta[jna]
            mi <- m[jna]
            ui <- S.inv[jna, -jna] %*% xo
            r <- (yi - xo %*% beta[-jna])[1, 1]
            taui <- tau[jna]

            # linear equation Ax = cc
            cc <- (r * betai / sigma^2 + mi - ui) / taui
            A <- (betai %*% t(betai) / sigma^2 + S.inv[jna, jna]) / (taui %*% t(taui))
            diag(A) <- 1

            # showEqn(A, cc)
            mu_tilde <- solve(A, cc)

            B_inv <- solve(A)

            Z <- mvrnorm(n = 1, mu = mu_tilde, Sigma = B_inv)
            X.sim[i, jna] <- Z / taui
          }
        }
      }
    }

    if (scale) {
      X.sim <- t(t(X.sim) * std.w)
      X.sim <- t(t(X.sim) + mean.w)
      mean.w <- apply(X.sim, 2, scale_mean.w, row.w)
      X.sim <- t(t(X.sim) - mean.w)
      std.w <- apply(X.sim, 2, scale_std.w, row.w) * sqrt(n)
      X.sim <- t(t(X.sim) / std.w)
    }

    # Stochastic Approximation step & Maximisation step
    beta.old <- beta
    sigma.old <- sigma
    theta.old <- theta
    if (missingcols != 0) {
      mu.old <- mu
      Sigma.old <- Sigma
    }

    # MLE of completed data likelihood
    # 1) beta
    W <- gamma * c + (rep(1, p) - gamma)
    revW <- 1 / W
    Xtemp <- sweep(X.sim, 2, revW, "*") # or: X.sim %*% diag(revW)
    z <- slope_admm(Xtemp, y, rep(0, p), rep(0, p), lambda_seq = lambda_sigma, 1)$z
    beta <- revW * z
    cstop <- sqrt(sum((beta - beta.old)^2))

    # 2) sigma
    rk <- p - rank(abs(z), ties.method = "max") + 1
    RSS <- sum((y - X.sim %*% beta)^2)
    sum_lamwbeta <- sum(lambda[rk] * abs(z))
    if (is.na(sigma.known)) {
      sigma <- (sqrt(sum_lamwbeta^2 + 4 * n * RSS) + sum_lamwbeta) / (2 * n)
    }

    # lambda
    lambda_sigma <- lambda * sigma
    lambda_sigma_inv <- lambda / sigma

    # 3) theta
    theta <- (sum(gamma == 1) + a) / (a + b + p)

    # 4) mu et Sigma
    if (missingcols != 0) {
      mu <- apply(X.sim, 2, mean)
      # Sigma = var(X.sim)*(n-1)/n
      Sigma <- linshrink_cov(X.sim) # linear shrinkage
    }

    if (eta != 1) {
      beta <- beta.old + eta * (beta - beta.old)
      sigma <- sigma.old + eta * (sigma - sigma.old)
      theta <- theta.old + eta * (theta - theta.old)
      if (missingcols != 0) {
        mu <- mu.old + eta * (mu - mu.old)
        Sigma <- Sigma.old + eta * (Sigma - Sigma.old)
      }
    }
    seqbeta[, t + 1] <- beta
    seqsigma[, t + 1] <- sigma
    seqgamma[, t + 1] <- gamma
    seqc[, t + 1] <- c
    seqtheta[, t + 1] <- theta
  }
  gamma.avg <- rowMeans(seqgamma[, -(1:20)], na.rm = TRUE)
  # if(!beta.scale) {beta.new = beta * (gamma.avg>1/2)}
  beta.new <- beta * (gamma.avg > 1 / 2)
  if (missingcols == 0) {
    mu <- Sigma <- NULL
  }
  return(list(X.sim = X.sim, beta = beta, seqbeta = seqbeta, beta.new = beta.new, gamma = gamma, gamma.avg = gamma.avg, seqgamma = seqgamma, pi.binom = pi.binom, theta = theta, seqtheta = seqtheta, sigma = sigma, seqsigma = seqsigma, c = c, seqc = seqc, mu = mu, Sigma = Sigma))
}

# Baseline of ABSLOPE
baseline_ABSLOPE <- function(n = 100, p = 20, nspr = 6, p.miss = 0.1, mu = rep(0, p), Sigma = diag(p), signallevel = 3, nb.seed = 1, a = 1, b = 10, tol_em = 1e-3, maxit = 300, impute = "mean", mec = "MCAR", sigma = 1, beta.start = NA, sigma.known = NA, sigma.init = NA, print_iter = FALSE, scale = FALSE, method_na = "lineq") {
  set.seed(nb.seed)
  # print(nb.seed)
  data.list <- data.generation(n, p, nspr, p.miss, mu, Sigma, signallevel, mec = mec, sigma = sigma)
  X <- data.list$X.obs
  X.comp <- data.list$X
  y <- data.list$y
  beta <- data.list$beta

  lambda <- create_lambda_bhq(ncol(X), fdr = 0.10)
  # start_time <- Sys.time()
  list.SLOB <- ABSLOPE(X, y, lambda, a = a, b = b, beta.start = beta.start, maxit = maxit, print_iter = print_iter, tol_em = tol_em, impute = impute, sigma.known = sigma.known, sigma.init = sigma.init, scale = scale, method_na = method_na)
  # time <- Sys.time() - start_time
  pr <- power(beta)
  fdr <- fdp(beta, which(list.SLOB$gamma.avg > 1 / 2))
  bias_beta_new <- norm(list.SLOB$beta.new - beta, "2") / norm(beta, "2")
  bias_sigma <- abs(list.SLOB$sigma - sigma) / sigma
  MSE_new <- norm(X.comp %*% list.SLOB$beta.new - X.comp %*% beta, "2") / norm(X.comp %*% beta, "2")
  return(list(pr = pr, fdr = fdr, bias_beta = bias_beta_new, bias_sigma = bias_sigma, MSE = MSE_new))
}
