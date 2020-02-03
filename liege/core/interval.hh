#ifndef LIEGE_CORE_INTERVAL_LOADED
#define LIEGE_CORE_INTERVAL_LOADED

// =============================================================================
// A simple interval for readible ranges such as cut parameters.
// =============================================================================
namespace liege {
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
} // namespace liege

#endif
