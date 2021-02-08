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

#ifndef LAGER_CORE_EXCEPTION_LOADED
#define LAGER_CORE_EXCEPTION_LOADED

#include <exception>
#include <string>

namespace lager {
class exception : public std::exception {
public:
  exception(const std::string& msg, const std::string& type = "exception")
      : msg_{msg}, type_{type} {}

  virtual const char* what() const throw() { return msg_.c_str(); }
  virtual const char* type() const throw() { return type_.c_str(); }
  virtual ~exception() throw() {}

private:
  std::string msg_;
  std::string type_;
};
} // ns lager

#endif
