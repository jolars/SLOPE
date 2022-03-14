#pragma once

#include "binomial.h"
#include "family.h"
#include "gaussian.h"
#include "multinomial.h"
#include "poisson.h"
#include <RcppArmadillo.h>
#include <memory>

// helper to choose family
template<typename... Ts>
std::unique_ptr<Family>
setupFamily(const std::string& family_choice, Ts... args)
{
  if (family_choice == "binomial")
    return std::unique_ptr<Binomial>(new Binomial{ std::forward<Ts>(args)... });
  else if (family_choice == "poisson")
    return std::unique_ptr<Poisson>(new Poisson{ std::forward<Ts>(args)... });
  else if (family_choice == "multinomial")
    return std::unique_ptr<Multinomial>(
      new Multinomial{ std::forward<Ts>(args)... });
  else
    return std::unique_ptr<Gaussian>(new Gaussian{ std::forward<Ts>(args)... });
}
