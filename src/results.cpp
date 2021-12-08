#include "results.h"
#include <RcppArmadillo.h>

Results::Results() {}

Results::Results(arma::mat beta,
                 arma::uword passes,
                 std::vector<double> primals,
                 std::vector<double> duals,
                 std::vector<double> time,
                 double deviance)
  : beta(beta)
  , passes(passes)
  , primals(primals)
  , duals(duals)
  , time(time)
  , deviance(deviance)
{}
