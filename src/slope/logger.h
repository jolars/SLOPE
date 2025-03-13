/**
 * @file
 * @brief Thread-safe warning logging facility for the slope library
 *
 * This module provides a centralized warning logging system that works
 * correctly in both single-threaded and multi-threaded (OpenMP) contexts.
 * Warnings are categorized by type and stored per thread to avoid contention.
 */

#pragma once

#include <map>
#include <mutex>
#include <string>

namespace slope {

/**
 * @brief Standard warning codes used throughout the slope library
 *
 * These enum values represent the various types of warnings that
 * can be generated during computational operations.
 */
enum class WarningCode
{
  GENERIC_WARNING,    ///< General uncategorized warning
  DEPRECATED_FEATURE, ///< Feature marked for deprecation
  MAXIT_REACHED,      ///< Maximum iterations reached without convergence
  LINE_SEARCH_FAILED  ///< Line search algorithm failed to converge
};

/**
 * @brief Convert a warning code to its string representation
 *
 * @param code The warning code to convert
 * @return std::string String representation of the warning code
 */
std::string
warningCodeToString(WarningCode code);

/**
 * @brief Thread-safe warning logger for tracking runtime warnings
 *
 * WarningLogger provides facilities for logging, retrieving, and managing
 * warnings that occur during computation. It is thread-safe and designed
 * to work with OpenMP parallel regions.
 *
 * Warnings are tracked per thread and can be accessed either individually
 * by thread or aggregated across all threads.
 */
class WarningLogger
{
private:
  /// Collection of warnings, organized by thread ID and warning code
  static std::map<int, std::map<WarningCode, std::string>> warnings;

  /// Mutex to protect concurrent access to the warnings collection
  static std::mutex warnings_mutex;

public:
  /**
   * @brief Log a new warning
   *
   * Stores a warning message associated with a specific warning code.
   * Thread-safe and can be called from within OpenMP parallel regions.
   * If a warning with the same code already exists for the current thread,
   * the message will be replaced.
   *
   * @param code The warning code categorizing this warning
   * @param message Descriptive message for the warning
   */
  static void addWarning(WarningCode code, const std::string& message);

  /**
   * @brief Retrieve all warnings from all threads
   *
   * Combines warnings from all threads into a single map. If multiple threads
   * logged warnings with the same code, the latest warning will be present
   * in the result.
   *
   * @return std::map<WarningCode, std::string> Map of all warning codes and
   * messages
   */
  static std::map<WarningCode, std::string> getWarnings();

  /**
   * @brief Clear all warnings
   *
   * Removes all warnings across all threads.
   */
  static void clearWarnings();

  /**
   * @brief Check if any warnings have been logged
   *
   * @return bool True if any thread has logged a warning, false otherwise
   */
  static bool hasWarnings();

  /**
   * @brief Get warnings from a specific thread
   *
   * Retrieves all warnings that were logged by a particular thread.
   *
   * @param threadId The ID of the thread whose warnings should be retrieved
   * @return std::map<WarningCode, std::string> Map of warning codes and
   * messages for the specified thread
   */
  static std::map<WarningCode, std::string> getThreadWarnings(int threadId);
};

} // namespace slope
