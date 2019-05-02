#ifndef PCSIM_CORE_LOGGER_LOADED
#define PCSIM_CORE_LOGGER_LOADED

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
namespace pcsim {
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
} // ns pcsim

// PREPROCESSOR macros to actually call the logger.
// Strongly prefered over calling the logger function directly, as in the macros,
// mtitle and mtext (which might be complex expressions) are only evaluated
// *after* the log_level check.
// This is *significantly* (orders of magnitude!) faster than calling 
// log<LEVEL>(mtitle, mtext) directly in the code.
#define LOG_CRITICAL(mtitle, mtext)                                            \
  if (pcsim::global::logger.level() >= log_level::CRITICAL) {               \
    pcsim::log<pcsim::log_level::CRITICAL>((mtitle), (mtext));           \
  }
#define LOG_ERROR(mtitle, mtext)                                               \
  if (pcsim::global::logger.level() >= log_level::ERROR) {                  \
    pcsim::log<pcsim::log_level::ERROR>((mtitle), (mtext));              \
  }
#define LOG_WARNING(mtitle, mtext)                                             \
  if (pcsim::global::logger.level() >= log_level::WARNING) {                \
    pcsim::log<pcsim::log_level::WARNING>((mtitle), (mtext));            \
  }
#define LOG_INFO(mtitle, mtext)                                                \
  if (pcsim::global::logger.level() >= log_level::INFO) {                   \
    pcsim::log<pcsim::log_level::INFO>((mtitle), (mtext));               \
  }
#define LOG_DEBUG(mtitle, mtext)                                               \
  if (pcsim::global::logger.level() >= log_level::DEBUG) {                  \
    pcsim::log<pcsim::log_level::DEBUG>((mtitle), (mtext));              \
  }
#define LOG_JUNK(mtitle, mtext)                                                \
  if (pcsim::global::logger.level() >= log_level::JUNK) {                   \
    pcsim::log<pcsim::log_level::JUNK>((mtitle), (mtext));               \
  }
#define LOG_JUNK2(mtitle, mtext)                                               \
  if (pcsim::global::logger.level() >= log_level::JUNK2) {                  \
    pcsim::log<pcsim::log_level::JUNK2>((mtitle), (mtext));              \
  }

// =============================================================================
// log_handler class designed for global usage,
// threading secure
// =============================================================================
namespace pcsim {
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
  void set_output(std::ostream& sink) { sink_ = &sink; }

  void operator()(const log_level mlevel, const std::string& mtitle,
                  const std::string& mtext) {
    lock_type lock{mutex_};
    if (mlevel > level_) {
      return;
    }
    time_t rt;
    time(&rt);

    std::ostream* sink =
        (mlevel == log_level::ERROR || mlevel == log_level::CRITICAL)
            ? &std::cerr
            : sink_;

    (*sink) << "[" << rt << ", " << mtitle << ", "
            << LOG_LEVEL_NAMES[static_cast<unsigned>(mlevel)] << "] " << mtext
            << std::endl;

    std::cout << "[" << rt << ", " << mtitle << ", "
            << LOG_LEVEL_NAMES[static_cast<unsigned>(mlevel)] << "] " << mtext
            << std::endl;
  }

private:
  log_level level_;
  std::ostream* sink_;
  mutable mutex_type mutex_;
};
} // ns pcsim

// =============================================================================
// global logger function
// threading seutil
// =============================================================================
namespace pcsim {
template <log_level level>
void log(const std::string& mtitle, const std::string& mtext,
         log_handler& logger) {
  logger(level, mtitle, mtext);
}
} // ns pcsim

#endif
