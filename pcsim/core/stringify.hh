#ifndef PCSIM_CORE_STRINGIFY_LOADED
#define PCSIM_CORE_STRINGIFY_LOADED

#include <string>
#include <utility>
#include <type_traits>

#include <boost/lexical_cast.hpp>

#include <pcsim/core/type_traits.hh>

namespace pcsim {

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
} // ns pcsim

// =============================================================================
// stringify(range, delimiter, accessor):
// return a stringified version of a range, with the elements separated by
// <delimiter>. The elements are accessed by <acc> prior to string conversion.
// Default accessor just returns the string, with specializations for STL
// containers and std::pair, allowing for nesting.
// =============================================================================
namespace pcsim {
template <class Range,
          class ElementAccessor = decltype(
              stringify_impl::element_accessor<typename Range::value_type>)>
std::string
stringify(const Range& r, const std::string& delimiter = ", ",
          ElementAccessor acc =
              stringify_impl::element_accessor<typename Range::value_type>);
} // ns pcsim

// =============================================================================
// Implementation
// =============================================================================
namespace pcsim {
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
} // ns pcsim

#endif
