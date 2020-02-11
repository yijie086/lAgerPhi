// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2020 Sylvester Joosten <sjoosten@anl.gov>
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

#ifndef LAGER_CORE_ASSERT_LOADED
#define LAGER_CORE_ASSERT_LOADED

#include <lager/core/logger.hh>
#include <lager/core/exception.hh>
#include <string>

// =============================================================================
// throwing assert throws an lager::exception with verbose error message.
// =============================================================================
#define tassert(condition, msg)                                                \
  if (!(condition)) {                                                          \
    lager::tassert_impl(#condition, __FILE__, __LINE__, msg);               \
  }

// =============================================================================
// implementation: throwing assert
// =============================================================================
namespace lager {
inline void tassert_impl(const std::string& condition,
                         const std::string& location, const int line,
                         const std::string& msg) {
  LOG_ERROR(location,
            "l" + std::to_string(line) + ": assert(" + condition + ") failed");
  throw lager::exception{msg, "assert"};
}
} // ns lager

#endif
