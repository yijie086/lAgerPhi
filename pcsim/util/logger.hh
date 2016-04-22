#ifndef PYTHIA6M_UTIL_LOGGER_LOADED
#define PYTHIA6M_UTIL_LOGGER_LOADED

#include <ctime>
#include <ostream>
#include <iostream> 
#include <string>
#include <array>
#include <mutex>

// =============================================================================
// default global logger
// =============================================================================

// =============================================================================
// Global logger, logging levels and functions.
//
// Actual logging should happen through the preprocessor macros
//   * LOG_CRITICAL(title, text)
//   * LOG_ERROR(title, text)
//   * LOG_WARNING(title, text)
//   * LOG_INFO(title, text)
//   * LOG_DEBUG(title, text)
//   * LOG_DEBUG2(title, text)
// =============================================================================
namespace pythia6m {
enum class log_level : unsigned {
  NOTHING = 0,
  CRITICAL = 1,
  ERROR = 2,
  WARNING = 3,
  INFO = 4,
  DEBUG = 5,
  JUNK = 6,
  JUNK2 = 7
};
constexpr std::array<const char*, 8> LOG_LEVEL_NAMES{
    "nothing", "critical", "error", "warning",
    "info",    "debug",    "junk",  "junk2"};

// the global logger
class log_handler;
namespace global {
extern log_handler logger;
} // ns global
template <log_level level>
void log(const std::string& mtitle, const std::string& mtext,
         log_handler& logger = global::logger);
} // ns pythia6m

// PREPROCESSOR macros to actually call the logger.
// Strongly prefered over calling the logger function directly, as in the macros,
// mtitle and mtext (which might be complex expressions) are only evaluated
// *after* the log_level check.
// This is *significantly* (orders of magnitude!) faster than calling 
// log<LEVEL>(mtitle, mtext) directly in the code.
#define LOG_CRITICAL(mtitle, mtext)                                            \
  if (pythia6m::global::logger.level() >= log_level::CRITICAL) {               \
    pythia6m::log<pythia6m::log_level::CRITICAL>((mtitle), (mtext));           \
  }
#define LOG_ERROR(mtitle, mtext)                                               \
  if (pythia6m::global::logger.level() >= log_level::ERROR) {                  \
    pythia6m::log<pythia6m::log_level::ERROR>((mtitle), (mtext));              \
  }
#define LOG_WARNING(mtitle, mtext)                                             \
  if (pythia6m::global::logger.level() >= log_level::WARNING) {                \
    pythia6m::log<pythia6m::log_level::WARNING>((mtitle), (mtext));            \
  }
#define LOG_INFO(mtitle, mtext)                                                \
  if (pythia6m::global::logger.level() >= log_level::INFO) {                   \
    pythia6m::log<pythia6m::log_level::INFO>((mtitle), (mtext));               \
  }
#define LOG_DEBUG(mtitle, mtext)                                               \
  if (pythia6m::global::logger.level() >= log_level::DEBUG) {                  \
    pythia6m::log<pythia6m::log_level::DEBUG>((mtitle), (mtext));              \
  }
#define LOG_JUNK(mtitle, mtext)                                                \
  if (pythia6m::global::logger.level() >= log_level::JUNK) {                   \
    pythia6m::log<pythia6m::log_level::JUNK>((mtitle), (mtext));               \
  }
#define LOG_JUNK2(mtitle, mtext)                                               \
  if (pythia6m::global::logger.level() >= log_level::JUNK2) {                  \
    pythia6m::log<pythia6m::log_level::JUNK2>((mtitle), (mtext));              \
  }

// =============================================================================
// log_handler class designed for global usage,
// threading seutil
// =============================================================================
namespace pythia6m {
class log_handler {
private:
  typedef std::mutex mutex_type;
  typedef std::lock_guard<mutex_type> lock_type;

public:
  log_handler(const log_level level = log_level::INFO,
              std::ostream& sink = std::cout);

  log_level level() const {
    lock_type lock{mutex_};
    return level_;
  }

  void set_level(const log_level level);
  void set_level(unsigned ulevel);

  void operator()(const log_level mlevel, const std::string& mtitle,
                  const std::string& mtext) {
    lock_type lock{mutex_};
    if (mlevel > level_)
      return;
    time_t rt;
    time(&rt);
    sink_ << "[" << rt << ", " << mtitle << ", "
          << LOG_LEVEL_NAMES[static_cast<unsigned>(mlevel)] << "] " << mtext
          << std::endl;
  }

private:
  log_level level_;
  std::ostream& sink_;
  mutable mutex_type mutex_;
};
} // ns pythia6m

// =============================================================================
// global logger function
// threading seutil
// =============================================================================
namespace pythia6m {
template <log_level level>
void log(const std::string& mtitle, const std::string& mtext,
         log_handler& logger) {
  logger(level, mtitle, mtext);
}
} // ns pythia6m

#endif
