/**
 * @file
 * @brief Simple high-resolution timer class for performance measurements
 */

#pragma once

#include <chrono>

namespace slope {

/**
 * @brief Timer class for measuring elapsed time with high resolution
 *
 * Uses std::chrono::high_resolution_clock to provide precise timing
 * measurements.
 */
class Timer
{
private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
  std::chrono::microseconds accumulated_time;
  bool is_running = false;

public:
  /**
   * @brief Starts the timer by recording the current time point
   */
  void start();

  /**
   * @brief Pauses the timer
   */
  void pause();

  /**
   * @brief Resumes the timer after a pause
   */
  void resume();

  /**
   * @brief Returns the elapsed time in seconds since start() was called
   * @return double Time elapsed in seconds with high precision
   */
  double elapsed() const;
};

} // namespace slope
