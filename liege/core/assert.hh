#ifndef LIEGE_CORE_ASSERT_LOADED
#define LIEGE_CORE_ASSERT_LOADED

#include <liege/core/logger.hh>
#include <liege/core/exception.hh>
#include <string>

// =============================================================================
// throwing assert throws an liege::exception with verbose error message.
// =============================================================================
#define tassert(condition, msg)                                                \
  if (!(condition)) {                                                          \
    liege::tassert_impl(#condition, __FILE__, __LINE__, msg);               \
  }

// =============================================================================
// implementation: throwing assert
// =============================================================================
namespace liege {
inline void tassert_impl(const std::string& condition,
                         const std::string& location, const int line,
                         const std::string& msg) {
  LOG_ERROR(location,
            "l" + std::to_string(line) + ": assert(" + condition + ") failed");
  throw liege::exception{msg, "assert"};
}
} // ns liege

#endif
