#include "utils.h"
#include <stdexcept>

namespace slope {

void
validateOption(const std::string& value,
               const std::set<std::string>& valid_options,
               const std::string& parameter_name)
{
  if (valid_options.find(value) == valid_options.end()) {
    std::string valid_list =
      std::accumulate(std::next(valid_options.begin()),
                      valid_options.end(),
                      std::string("'") + *valid_options.begin() + "'",
                      [](const std::string& a, const std::string& b) {
                        return a + ", '" + b + "'";
                      });

    throw std::invalid_argument("Invalid " + parameter_name + ": '" + value +
                                "'. Must be one of: " + valid_list);
  }
}

} // namespace slope
