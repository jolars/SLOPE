firstUpper <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

camelCase <- function(x) {
  s <- strsplit(x, "[^[:alnum:]]")

  sapply(s, function(y) {
    first <- toupper(substring(y, 1, 1))
    paste(first, substring(y, 2), sep = "", collapse = "")
  })
}

randomProblem <- function(n = 1000,
                          p = 100,
                          q = 0.2,
                          n_groups = NULL,
                          n_targets =
                            if (match.arg(response) == "multinomial") 3 else 1,
                          density = 1,
                          amplitude =
                            if (match.arg(response) == "poisson") 1 else 3,
                          alpha = 1,
                          response = c(
                            "gaussian",
                            "binomial",
                            "poisson",
                            "multinomial"
                          ),
                          rho = 0) {
  m <- n_targets

  if (density == 1) {
    x <- matrix(stats::rnorm(n * p), n)
  } else {
    x <- Matrix::rsparsematrix(n, p, density)
  }

  if (rho > 0) {
    x <- sqrt(1 - rho) * x + sqrt(rho) * stats::rnorm(n)
  }

  if (!is.null(n_groups)) {
    groups <- rep(seq_len(n_groups),
      each = ceiling(m * p / n_groups),
      length.out = p * m
    )
    nonzero <- which(groups %in% seq_len(max(floor(n_groups * q), 1)))
  } else {
    groups <- NA
    nonzero <- sample(p * m, max(floor(q * p * m), 1))
  }

  signs <- sample(c(-1, 1), p * m, replace = TRUE)

  beta <- signs * amplitude * (1:(p * m) %in% nonzero)

  y <- switch(
    match.arg(response),
    gaussian = x %*% beta + stats::rnorm(n, sd = alpha),
    binomial = {
      y <- x %*% beta + stats::rnorm(n, sd = alpha)
      (sign(y) + 1) / 2
    },
    multinomial = {
      beta <- matrix(beta, p, m)
      lin_pred_exp <- exp(x %*% beta)
      prob <- lin_pred_exp / rowSums(lin_pred_exp)
      y <- apply(prob, 1, function(x) sample(1:m, 1, prob = x))
    },
    poisson = {
      lambda <- as.double(exp(x %*% beta))
      y <- stats::rpois(n, lambda)
    }
  )

  dimnames(x) <- list(
    seq_len(nrow(x)),
    paste0("V", seq_len(ncol(x)))
  )

  list(
    x = x,
    y = as.double(y),
    beta = beta,
    groups = groups,
    nonzero = nonzero,
    q = q
  )
}
