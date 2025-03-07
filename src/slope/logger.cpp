#include "logger.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace slope {

// Define static member variables
std::map<int, std::map<WarningCode, std::string>> WarningLogger::warnings;
std::mutex WarningLogger::warnings_mutex;

std::string
warningCodeToString(WarningCode code)
{
  switch (code) {
    case WarningCode::GENERIC_WARNING:
      return "GENERIC_WARNING";
    case WarningCode::DEPRECATED_FEATURE:
      return "DEPRECATED_FEATURE";
    case WarningCode::MAXIT_REACHED:
      return "MAXIT_REACHED";
    case WarningCode::LINE_SEARCH_FAILED:
      return "LINE_SEARCH_FAILED";
    default:
      return "UNKNOWN_WARNING";
  }
}

void
WarningLogger::addWarning(WarningCode code, const std::string& message)
{
  int threadId = 0;
#ifdef _OPENMP
  threadId = omp_get_thread_num();
#endif

  std::lock_guard<std::mutex> lock(warnings_mutex);
  warnings[threadId][code] = message;
}

std::map<WarningCode, std::string>
WarningLogger::getWarnings()
{
  std::lock_guard<std::mutex> lock(warnings_mutex);
  std::map<WarningCode, std::string> allWarnings;
  for (const auto& threadWarnings : warnings) {
    allWarnings.insert(threadWarnings.second.begin(),
                       threadWarnings.second.end());
  }
  return allWarnings;
}

void
WarningLogger::clearWarnings()
{
  std::lock_guard<std::mutex> lock(warnings_mutex);
  warnings.clear();
}

bool
WarningLogger::hasWarnings()
{
  std::lock_guard<std::mutex> lock(warnings_mutex);
  for (const auto& threadWarnings : warnings) {
    if (!threadWarnings.second.empty()) {
      return true;
    }
  }
  return false;
}

std::map<WarningCode, std::string>
WarningLogger::getThreadWarnings(int threadId)
{
  std::lock_guard<std::mutex> lock(warnings_mutex);
  if (warnings.find(threadId) != warnings.end()) {
    return warnings[threadId];
  }
  return std::map<WarningCode, std::string>();
}

} // namespace slope
