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

#ifndef LAGER_CORE_STRINGIFY_LOADED
#define LAGER_CORE_STRINGIFY_LOADED

#include <boost/lexical_cast.hpp>
#include <lager/core/type_traits.hh>
#include <string>
#include <type_traits>
#include <utility>

namespace lager {

// implementation prototypes needed for stringify
namespace stringify_impl {
// standard accessor returns el untouched
template <class RangeElement, class = typename std::enable_if<
                                  !is_container<RangeElement>::value>::type>
RangeElement element_accessor(const RangeElement& el);
// container accessor returns stringified version of container
template <class RangeElement, class = typename std::enable_if<
                                  is_container<RangeElement>::value>::type,
          class = void>
RangeElement element_accessor(const RangeElement& el);
// pair element accessor returns a "<first>: <second>" string
template <class U, class V>
std::string element_accessor(const std::pair<U, V>& p);
}
} // ns lager

// =============================================================================
// stringify(range, delimiter, accessor):
// return a stringified version of a range, with the elements separated by
// <delimiter>. The elements are accessed by <acc> prior to string conversion.
// Default accessor just returns the string, with specializations for STL
// containers and std::pair, allowing for nesting.
// =============================================================================
namespace lager {
template <class Range,
          class ElementAccessor = decltype(
              stringify_impl::element_accessor<typename Range::value_type>)>
std::string
stringify(const Range& r, const std::string& delimiter = ", ",
          ElementAccessor acc =
              stringify_impl::element_accessor<typename Range::value_type>);
} // ns lager

// =============================================================================
// Implementation
// =============================================================================
namespace lager {
namespace stringify_impl {
template <class RangeElement, class>
RangeElement element_accessor(const RangeElement& el) {
  return el;
}
template <class RangeElement, class, class>
RangeElement element_accessor(const RangeElement& el) {
  return "(" + stringify(el) + ")";
}
template <class U, class V>
std::string element_accessor(const std::pair<U, V>& p) {
  return (boost::lexical_cast<std::string>(element_accessor(p.first)) + ": " +
          boost::lexical_cast<std::string>(element_accessor(p.second)));
}
}

template <class Range, class ElementAccessor>
std::string stringify(const Range& r, const std::string& delimiter,
                      ElementAccessor acc) {
  std::string str{};
  for (const auto& el : r) {
    if (!str.empty()) {
      str += delimiter;
    }
    str += boost::lexical_cast<std::string>(acc(el));
  }
  return str;
}
} // ns lager

#endif
