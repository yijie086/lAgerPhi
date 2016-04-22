#ifndef PYTHIA6M_UTIL_ASSERT_LOADED
#define PYTHIA6M_UTIL_ASSERT_LOADED

#include <pythia6m/util/logger.hh>
#include <pythia6m/util/exception.hh>
#include <string>

// =============================================================================
// throwing assert throws an pythia6m::exception with verbose error message.
// =============================================================================
#define tassert(condition, msg)                                                \
  if (!(condition)) {                                                          \
    pythia6m::tassert_impl(#condition, __FILE__, __LINE__, msg);               \
  }

// =============================================================================
// implementation: throwing assert
// =============================================================================
namespace pythia6m {
inline void tassert_impl(const std::string& condition,
                         const std::string& location, const int line,
                         const std::string& msg) {
  LOG_ERROR(location,
            "l" + std::to_string(line) + ": assert(" + condition + ") failed");
  throw pythia6m::exception{msg, "assert"};
}
} // ns pythia6m

#endif
