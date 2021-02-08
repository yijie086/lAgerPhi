// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
// 
// This file is part of lAger.
// 
// lAger is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Shoftware Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// lAger is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with lAger.  If not, see <https://www.gnu.org/licenses/>.
// 

#ifndef LAGER_CORE_LOGGER_LOADED
#define LAGER_CORE_LOGGER_LOADED

#include <array>
#include <ctime>
#include <iostream>
#include <mutex>
#include <ostream>
#include <string>

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
namespace lager {
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
} // ns lager

// PREPROCESSOR macros to actually call the logger.
// Strongly prefered over calling the logger function directly, as in the
// macros,
// mtitle and mtext (which might be complex expressions) are only evaluated
// *after* the log_level check.
// This is *significantly* (orders of magnitude!) faster than calling
// log<LEVEL>(mtitle, mtext) directly in the code.
#define LOG_CRITICAL(mtitle, mtext)                                            \
  if (lager::global::logger.level() >= log_level::CRITICAL) {                  \
    lager::log<lager::log_level::CRITICAL>((mtitle), (mtext));                 \
  }
#define LOG_ERROR(mtitle, mtext)                                               \
  if (lager::global::logger.level() >= log_level::ERROR) {                     \
    lager::log<lager::log_level::ERROR>((mtitle), (mtext));                    \
  }
#define LOG_WARNING(mtitle, mtext)                                             \
  if (lager::global::logger.level() >= log_level::WARNING) {                   \
    lager::log<lager::log_level::WARNING>((mtitle), (mtext));                  \
  }
#define LOG_INFO(mtitle, mtext)                                                \
  if (lager::global::logger.level() >= log_level::INFO) {                      \
    lager::log<lager::log_level::INFO>((mtitle), (mtext));                     \
  }
#define LOG_DEBUG(mtitle, mtext)                                               \
  if (lager::global::logger.level() >= log_level::DEBUG) {                     \
    lager::log<lager::log_level::DEBUG>((mtitle), (mtext));                    \
  }
#define LOG_JUNK(mtitle, mtext)                                                \
  if (lager::global::logger.level() >= log_level::JUNK) {                      \
    lager::log<lager::log_level::JUNK>((mtitle), (mtext));                     \
  }
#define LOG_JUNK2(mtitle, mtext)                                               \
  if (lager::global::logger.level() >= log_level::JUNK2) {                     \
    lager::log<lager::log_level::JUNK2>((mtitle), (mtext));                    \
  }

// =============================================================================
// log_handler class designed for global usage,
// threading secure
// =============================================================================
namespace lager {
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

    if (sink != &std::cout) {
      (*sink) << "[" << rt << ", " << mtitle << ", "
              << LOG_LEVEL_NAMES[static_cast<unsigned>(mlevel)] << "] " << mtext
              << std::endl;
    }

    std::cout << "[" << rt << ", " << mtitle << ", "
              << LOG_LEVEL_NAMES[static_cast<unsigned>(mlevel)] << "] " << mtext
              << std::endl;
  }

private:
  log_level level_;
  std::ostream* sink_;
  mutable mutex_type mutex_;
};
} // ns lager

// =============================================================================
// global logger function
// threading seutil
// =============================================================================
namespace lager {
template <log_level level>
void log(const std::string& mtitle, const std::string& mtext,
         log_handler& logger) {
  logger(level, mtitle, mtext);
}
} // ns lager

#endif
