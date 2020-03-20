#' Setup a data.frame of diagnostics
#'
#' @param res the result from calling the C++ routine used to fit a model
#'   in SLOPE
#'
#' @return A data.frame
#'
#' @keywords internal
setupDiagnostics <- function(res) {
  time <- res$time
  primals <- res$primals
  duals <- res$duals

  nl <- length(time)
  nn <- lengths(time)
  time <- unlist(time)
  primal <- unlist(primals)
  dual <- unlist(duals)

  data.frame(iteration = unlist(lapply(nn, seq_len)),
             time = time,
             primal = primal,
             dual = dual,
             penalty = rep(seq_len(nl), nn))
}

