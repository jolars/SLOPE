#include "timer.h"

namespace slope {

void
Timer::start()
{
  start_time = std::chrono::high_resolution_clock::now();
}

double
Timer::elapsed() const
{
  auto current_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
    current_time - start_time);

  return duration.count() * 1e-6; // Convert microseconds to seconds
}

} // namespace slope
