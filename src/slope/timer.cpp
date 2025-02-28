#include "timer.h"

namespace slope {

void
Timer::start()
{
  start_time = std::chrono::high_resolution_clock::now();
  is_running = true;
  accumulated_time = std::chrono::microseconds(0);
}

void
Timer::pause()
{
  if (is_running) {
    auto current_time = std::chrono::high_resolution_clock::now();
    accumulated_time += std::chrono::duration_cast<std::chrono::microseconds>(
      current_time - start_time);
    is_running = false;
  }
}

void
Timer::resume()
{
  if (!is_running) {
    start_time = std::chrono::high_resolution_clock::now();
    is_running = true;
  }
}

double
Timer::elapsed() const
{
  auto total_time = accumulated_time;
  if (is_running) {
    auto current_time = std::chrono::high_resolution_clock::now();
    total_time += std::chrono::duration_cast<std::chrono::microseconds>(
      current_time - start_time);
  }
  return total_time.count() * 1e-6; // Convert microseconds to seconds
}

} // namespace slope
