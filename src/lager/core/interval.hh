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

#ifndef LAGER_CORE_INTERVAL_LOADED
#define LAGER_CORE_INTERVAL_LOADED

// =============================================================================
// A simple interval for readible ranges such as cut parameters.
// =============================================================================
namespace lager {
template <class T> struct interval {
  using value_type = T;
  value_type min;
  value_type max;
  constexpr interval() = default;
  constexpr interval(const value_type min, const value_type max)
      : min{min}, max{max} {}
  constexpr interval(const std::pair<value_type, value_type>& rhs)
      : min{rhs.first}, max{rhs.second} {}
  constexpr bool includes(const value_type value) const {
    return value > min && value < max;
  }
  constexpr bool excludes(const value_type value) const {
    return !includes(value);
  }
  constexpr T width() const { return max - min; }
  // Implicitly conversion to std::pair<T,T> when needed.
  constexpr operator std::pair<value_type, value_type>() const {
    return std::make_pair(min, max);
  }
  constexpr interval operator*(const double rhs) const {
    return {min * rhs, max * rhs};
  }
  constexpr interval operator/(const double rhs) const {
    return operator*(1 / rhs);
  }
  interval operator*=(const double rhs) { *this = operator*(rhs); }
  interval operator/=(const double rhs) { *this = operator/(rhs); }
};
} // namespace lager

#endif
