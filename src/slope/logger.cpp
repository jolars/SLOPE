#include "logger.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace slope {

// Define static member variables
std::map<int, std::vector<Warning>> WarningLogger::warnings;
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
  int thread_id = 0;
#ifdef _OPENMP
  thread_id = omp_get_thread_num();
#endif

  std::lock_guard<std::mutex> lock(warnings_mutex);
  warnings[thread_id].emplace_back(code, message);
}

std::vector<Warning>
WarningLogger::getWarnings()
{
  std::lock_guard<std::mutex> lock(warnings_mutex);
  std::vector<Warning> all_warnings;
  for (const auto& thread_warnings : warnings) {
    all_warnings.insert(all_warnings.end(),
                        thread_warnings.second.begin(),
                        thread_warnings.second.end());
  }
  return all_warnings;
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
  for (const auto& thread_warnings : warnings) {
    if (!thread_warnings.second.empty()) {
      return true;
    }
  }

  return false;
}
std::vector<Warning>
WarningLogger::getThreadWarnings(int thread_id)
{
  std::lock_guard<std::mutex> lock(warnings_mutex);
  if (warnings.find(thread_id) != warnings.end()) {
    return warnings[thread_id];
  }
  return std::vector<Warning>();
}

} // namespace slope
