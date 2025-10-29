#include <slope/losses/logistic.h>
#include <slope/losses/multinomial.h>
#include <slope/losses/poisson.h>
#include <slope/losses/quadratic.h>
#include <slope/losses/setup_loss.h>

namespace slope {

std::unique_ptr<Loss>
setupLoss(const std::string& loss)
{
  if (loss == "logistic")
    return std::make_unique<Logistic>();
  else if (loss == "poisson")
    return std::make_unique<Poisson>();
  else if (loss == "multinomial")
    return std::make_unique<Multinomial>();

  // else Quadratic
  return std::make_unique<Quadratic>();
}

}
