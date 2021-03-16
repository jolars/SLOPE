preprocessResponse <- function(family, y) {
  switch(
    family,
    gaussian = {
      y <- as.numeric(y)

      if (NCOL(y) > 1)
        stop("response for Gaussian regression must be one-dimensional.")

      y_center <- mean(y)
      y_scale  <- 1

      y <- as.matrix(y - y_center)

      list(y = y,
           y_center = y_center,
           y_scale = y_scale,
           n_classes = 1L,
           n_targets = 1L,
           class_names = NA_character_)
    },

    binomial = {
      if (NCOL(y) > 1)
        stop("response for binomial regression must be one-dimensional.")

      if (length(unique(y)) > 2)
        stop("more than two classes in response")

      if (length(unique(y)) == 1)
        stop("only one class in response.")

      y_table <- table(y)
      min_class <- min(y_table)

      if (min_class <= 1)
        stop("one class only has ", min_class, " observations.")

      class_names <- names(y_table)

      # Transform response to {-1, 1}, which is used internally
      y <- as.matrix(ifelse(as.numeric(as.factor(y)) == 1, -1, 1))

      list(y = y,
           y_center = 0,
           y_scale = 1,
           n_classes = 1L,
           n_targets = 1L,
           class_names = class_names,
           response_names = colnames(y))
    },

    multinomial = {
      if (NCOL(y) > 1)
        stop("response for multinomial regression must be one-dimensional.")

      y <- droplevels(as.factor(y))
      y_table <- table(y)
      min_class <- min(y_table)
      class_names <- names(y_table)
      n_classes <- length(y_table)
      n_targets <- n_classes - 1

      Y <- matrix(NA, NROW(y), n_targets)

      for (k in 1:n_targets) {
        Y[, k] <- as.integer(y == names(y_table)[k])
      }

      if (n_classes == 2)
        stop("only two classes in response. Are you looking for family = 'binomial'?")

      if (n_classes == 1)
        stop("only one class in response")

      if (min_class <= 1)
        stop("one class only has ", min_class, " observations.")

      list(y = Y,
           y_center = rep(0, n_targets),
           y_scale = rep(1, n_targets),
           n_classes = n_classes,
           n_targets = n_targets,
           class_names = class_names,
           response_names = class_names)
    },

    poisson = {
      if (NCOL(y) > 1)
        stop("response for poisson regression must be one-dimensional.")

      if (any(y < 0))
        stop("cannot have negative responses in poisson model")

      list(y = y,
           y_center = 0,
           y_scale = 1,
           n_classes = 1L,
           n_targets = 1L,
           class_names = NA_character_,
           response_names = colnames(y))
    }
  )
}
