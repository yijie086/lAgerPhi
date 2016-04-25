#ifndef PCSIM_UTIL_ASSERT_LOADED
#define PCSIM_UTIL_ASSERT_LOADED

#include <pcsim/util/logger.hh>
#include <pcsim/util/exception.hh>
#include <string>

// =============================================================================
// throwing assert throws an pcsim::exception with verbose error message.
// =============================================================================
#define tassert(condition, msg)                                                \
  if (!(condition)) {                                                          \
    pcsim::tassert_impl(#condition, __FILE__, __LINE__, msg);               \
  }

// =============================================================================
// implementation: throwing assert
// =============================================================================
namespace pcsim {
inline void tassert_impl(const std::string& condition,
                         const std::string& location, const int line,
                         const std::string& msg) {
  LOG_ERROR(location,
            "l" + std::to_string(line) + ": assert(" + condition + ") failed");
  throw pcsim::exception{msg, "assert"};
}
} // ns pcsim

#endif
