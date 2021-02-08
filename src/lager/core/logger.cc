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

#include "logger.hh"

namespace lager {
// =============================================================================
// Implementation: log_handler
// =============================================================================
log_handler::log_handler(const log_level level, std::ostream& sink)
    : level_{level}, sink_(&sink) {}

void log_handler::set_level(const log_level level) {
  lock_type lock{mutex_};
  level_ = level;
}

void log_handler::set_level(unsigned ulevel) {
  lock_type lock{mutex_};
  if (ulevel >= LOG_LEVEL_NAMES.size()) {
    ulevel = LOG_LEVEL_NAMES.size() - 1;
  }
  level_ = static_cast<log_level>(ulevel);
}

// =============================================================================
// instantiation of the global logger 
// =============================================================================
namespace global {
log_handler logger{};
} // ns global
} // ns lager
