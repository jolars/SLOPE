#' @keywords internal
#' @aliases SLOPE-package
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @useDynLib SLOPE, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import Matrix
#' @import foreach
## usethis namespace: end

if(getRversion() >= "2.15.1")
  utils::globalVariables(c('Freq', 'Var1', 'alpha', 'lo', 'hi'))
