#ifndef PCSIM_CORE_CONFIGURATION_LOADED
#define PCSIM_CORE_CONFIGURATION_LOADED

#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/property_tree/exceptions.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <map>
#include <pcsim/core/exception.hh>
#include <pcsim/core/interval.hh>
#include <pcsim/core/stringify.hh>
#include <string>

namespace pcsim {
// necessary type aliases
using ptree = boost::property_tree::ptree;
using boost::property_tree::read_json;
using boost::property_tree::write_json;
template <class T> using optional = boost::optional<T>;
template <class T> using translation_map = std::map<std::string, T>;

// error prototypes
class configuration_error;
class configuration_path_error;
class configuration_key_error;
class configuration_value_error;
class configuration_translation_error;
} // ns pcsim

// =============================================================================
// Consistent way to format object paths/names/titles
//
// Note: * for optimal compatiblity with ROOT, a '/' path separator should be
//         used.
//       * a format_path() is not provided, as this is performed more
//         transparently by the string_path class.
// =============================================================================
namespace pcsim {
// path/name/title separators.
constexpr const char PATH_SEPARATOR{'/'};
constexpr const char NAME_SEPARATOR{'_'};
constexpr const char* const TITLE_SEPARATOR{" "};

// utility functions to consistently format object names and titles
std::string format_name(std::string context, const std::string& name);
std::string format_title(const std::string& context, std::string title);

// Utility class to store string paths (based on ptree::path_type),
// automatically taking care of the prefered PATH_SEPARATOR.
// Casts to std::string to make dumping the path contents less cumbersome.
class string_path : public ptree::path_type {
public:
  using base_type = ptree::path_type;
  explicit string_path(const char separator = PATH_SEPARATOR)
      : base_type{separator} {}
  string_path(const std::string& path, const char separator = PATH_SEPARATOR)
      : base_type{path, separator} {}
  string_path(const char* path, const char separator = PATH_SEPARATOR)
      : string_path(std::string{path}, separator) {}
  string_path(const base_type& rhs) : base_type{rhs} {}
  string_path& operator=(const string_path& rhs) {
    static_cast<base_type>(*this) = rhs;
    return *this;
  }
  std::string str() const { return dump(); }
};
} // ns pcsim

// =============================================================================
// configuration handler
//
// The configuration constructor takes an identifier string
// as argument. This string should match the settings path
// in the associated settings ptree.
// =============================================================================
namespace pcsim {

namespace configuration_impl {
  class value_proxy;
}

class configuration {
public:
  constexpr static const char* DEFAULTS{"defaults"};
  constexpr static const char* TYPE_KEY{"type"};

  configuration(const ptree& settings, const string_path& path);
  configuration(const configuration& conf, const string_path& path)
      : configuration{conf.settings_, path} {}

  // load the settings from a given ptree
  void load(const ptree& in_conf);

  // store the settings in the give ptree
  // The defaults are only exported if not yet present.
  void save(ptree& out_conf) const;

  // get the type info
  std::string type() const { return get<std::string>(TYPE_KEY); }

  // add a value to the configuration
  template <class T> void set(const string_path& key, const T& value) {
    settings_.put(key, value);
  }

  // Universal getter
  // Use these instead of the get<T>() member functions below for more readible
  // code
  configuration_impl::value_proxy value(const string_path& key) const;
  configuration_impl::value_proxy operator[](const string_path& key) const;

  // Three pairs of functions to get a setting by its key.
  //
  // In each pair, the translation_map version will lookup the
  // configuration value in the map, and throw a
  // configuration_translation_error if the lookup failed.
  //
  // 1. optional version
  template <class T> optional<T> get_optional(const string_path& key) const;
  template <class T>
  optional<T> get_optional(const string_path& key,
                           const translation_map<T>& tr) const;
  // 2. Throwing version
  template <class T> T get(const string_path& key) const;
  template <class T>
  T get(const string_path& key, const translation_map<T>& tr) const;
  // 3. default-value version
  //      if default_value is needed, it is automatically added to
  //      the default configurations for this entity
  template <class T, class = typename std::enable_if<!is_map<T>::value>::type>
  T get(const string_path& key, const T& default_value);
  template <class T>
  T get(const string_path& key, const T& default_value,
        const translation_map<T>& tr);
  // same getters, but in (STL) vector version
  // 1. optional version
  template <class T>
  optional<std::vector<T>> get_optional_vector(const string_path& key) const;
  template <class T>
  optional<std::vector<T>>
  get_optional_vector(const string_path& key,
                      const translation_map<T>& tr) const;
  // 2. Throwing version
  template <class T> std::vector<T> get_vector(const string_path& key) const;
  template <class T>
  std::vector<T> get_vector(const string_path& key,
                            const translation_map<T>& tr) const;

  // special physics vector versions
  // 1. optional versions
  template <class Vector3>
  optional<Vector3> get_optional_vector3(const string_path& key) const;
  template <class Vector3>
  optional<Vector3>
  get_optional_vector3(const string_path& key,
                       const translation_map<double>& tr) const;
  template <class Vector4>
  optional<Vector4> get_optional_vector4(const string_path& key) const;
  template <class Vector4>
  optional<Vector4>
  get_optional_vector4(const string_path& key,
                       const translation_map<double>& tr) const;
  // 2. throwing versions
  template <class Vector3> Vector3 get_vector3(const string_path& key) const;
  template <class Vector3>
  Vector3 get_vector3(const string_path& key,
                      const translation_map<double>& tr) const;
  template <class Vector4> Vector4 get_vector4(const string_path& key) const;
  template <class Vector4>
  Vector4 get_vector4(const string_path& key,
                      const translation_map<double>& tr) const;

  // special version to create bit pattern of a vector of
  // bit patterns
  // 1. optional versions
  template <class T>
  optional<T> get_optional_bitpattern(const string_path& key) const;
  template <class T>
  optional<T> get_optional_bitpattern(const string_path& key,
                                      const translation_map<T>& tr) const;
  // 2. throwing versions
  template <class T> T get_bitpattern(const string_path& key) const;
  template <class T>
  T get_bitpattern(const string_path& key, const translation_map<T>& tr) const;
  // special version to get an interval from a "range" vector
  // 1. optional versions
  template <class T>
  optional<interval<T>> get_optional_range(const string_path& key) const;
  template <class T>
  optional<interval<T>> get_optional_range(const string_path& key,
                                           const translation_map<T>& tr) const;
  // 2. throwing versions
  template <class T> interval<T> get_range(const string_path& key) const;
  template <class T>
  interval<T> get_range(const string_path& key,
                        const translation_map<T>& tr) const;


  // Helper functions to construct exceptions
  configuration_path_error path_error(const string_path& path) const;
  configuration_key_error key_error(const string_path& key) const;
  configuration_value_error value_error(const string_path& key,
                                        const std::string& value) const;
  configuration_value_error value_error(const string_path& key) const;
  configuration_translation_error
  translation_error(const string_path& key, const std::string& value) const;
  configuration_translation_error
  translation_error(const string_path& key) const;
  template <class T>
  configuration_translation_error
  translation_error(const string_path& key, const std::string& value,
                    const translation_map<T>& tr) const;

private:
  template <class T>
  T translate(const string_path& key, const std::string& val,
              const translation_map<T>& tr) const;

  // settings
  string_path settings_path_;
  ptree settings_;
  // defaults
  string_path defaults_path_;
  ptree defaults_;
};
} // ns pcsim

// =============================================================================
// configurable mixin
//
// Derive classes that should be configurable from the configurable mixin to
// automatically take care of configuration storage/loading/saving.
//
// The path member function returns the configuration path of the object, making
// it easier to create path hierarchies for nested objects.
// =============================================================================
namespace pcsim {
class configurable {
public:
  // ptree-constructor is deprecated, configuration constructor prefered
  configurable(const ptree& settings, const string_path& path)
      : conf_{settings, path}, path_{path} {}
  configurable(const configuration& conf, const string_path& path)
      : conf_{conf, path}, path_{path} {}
  const string_path& path() const { return path_; }
  const configuration& conf() const { return conf_; }
  configuration& conf() { return conf_; }

private:
  configuration conf_;
  const string_path path_;
};
} // ns pcsim


// =============================================================================
// Definition: exceptions
// =============================================================================
namespace pcsim {
// exceptions
class configuration_error : public pcsim::exception {
public:
  configuration_error(const std::string& msg,
                      const std::string& type = "configuration_error")
      : pcsim::exception{msg, type} {}
};

class configuration_path_error : public configuration_error {
public:
  configuration_path_error(const string_path& path);
};
class configuration_key_error : public configuration_error {
public:
  configuration_key_error(const string_path& key,
                          const string_path& settings_path,
                          const string_path& defaults_path);
};

class configuration_value_error : public configuration_error {
public:
  configuration_value_error(const string_path& key, const std::string& value,
                            const string_path& settings_path,
                            const string_path& defaults_path);
};
class configuration_translation_error : public configuration_error {
public:
  configuration_translation_error(const string_path& key,
                                  const std::string& value,
                                  const string_path& settings_path,
                                  const string_path& defaults_path);
  template <class T>
  configuration_translation_error(const string_path& key,
                                  const std::string& value,
                                  const translation_map<T>& tr,
                                  const string_path& settings_path,
                                  const string_path& defaults_path);
};

// =============================================================================
// Implementation: configuration value proxies with auto casting
// =============================================================================
namespace configuration_impl {
class value_proxy {
public:
  value_proxy(const configuration& conf, string_path key)
      : conf_{conf}, key_{std::move(key)} {}

  template <class T> operator T() const { return conf_.get<T>(key_); }
  template <class T> operator std::vector<T>() const {
    return conf_.get_vector<T>(key_);
  }
  template <class T> operator interval<T>() const {
    return conf_.get_range<T>(key_);
  }

private:
  const configuration& conf_;
  const string_path key_;
};
}

// =============================================================================
// Implementation: configuration getters
// =============================================================================
template <class T>
optional<T> configuration::get_optional(const string_path& key) const {
  try {
    auto s = settings_.get_optional<T>(key);
    if (!s) {
      s = defaults_.get_optional<T>(key);
    }
    return s;
  } catch (boost::property_tree::ptree_bad_data& e) {
    throw translation_error(key, e.data<std::string>());
  }
}
template <class T>
optional<T> configuration::get_optional(const string_path& key,
                                        const translation_map<T>& tr) const {
  auto s = get_optional<std::string>(key);
  if (!s) {
    return {};
  }
  return {translate(key, *s, tr)};
}
template <class T> T configuration::get(const string_path& key) const {
  auto s = get_optional<T>(key);
  if (!s) {
    throw key_error(key);
  }
  return *s;
}
template <class T>
T configuration::get(const string_path& key,
                     const translation_map<T>& tr) const {
  std::string val{get<std::string>(key)};
  return translate(key, val, tr);
}
template <class T, class>
T configuration::get(const string_path& key, const T& default_value) {
  auto s = get_optional<T>(key);
  if (!s) {
    defaults_.put(key, default_value);
    return default_value;
  }
  return *s;
}
template <class T>
T configuration::get(const string_path& key, const T& default_value,
                     const translation_map<T>& tr) {
  std::string val{get(key, default_value)};
  return translate(key, val, tr);
}
// and vector versions
template <class T>
optional<std::vector<T>>
configuration::get_optional_vector(const string_path& key) const {
  optional<std::vector<T>> vec;
  auto node = settings_.get_child_optional(key);
  if (!node) {
    node = defaults_.get_child_optional(key);
  }
  if (node) {
    vec.reset(std::vector<T>());
    for (const auto& child : *node) {
      auto val = child.second.get_value_optional<T>();
      if (val) {
        vec->push_back(*val);
      }
    }
  }
  return vec;
}
template <class T>
optional<std::vector<T>>
configuration::get_optional_vector(const string_path& key,
                                   const translation_map<T>& tr) const {
  optional<std::vector<T>> vec;
  auto vec_str = get_optional_vector<std::string>(key);
  if (vec_str) {
    vec.reset(std::vector<T>());
    for (const auto& el : *vec_str) {
      vec->push_back(translate(key, el, tr));
    }
  }
  return vec;
}
template <class T>
std::vector<T> configuration::get_vector(const string_path& key) const {
  auto s = get_optional_vector<T>(key);
  if (!s) {
    throw key_error(key);
  }
  return *s;
}
template <class T>
std::vector<T> configuration::get_vector(const string_path& key,
                                         const translation_map<T>& tr) const {
  auto s = get_optional_vector<T>(key, tr);
  if (!s) {
    throw key_error(key);
  }
  return *s;
}
// and 3-Vector version
template <class Vector3>
optional<Vector3>
configuration::get_optional_vector3(const string_path& key) const {
  auto range = get_optional_vector<double>(key);
  if (range) {
    if (range->size() != 3) {
      throw translation_error(key, stringify(*range));
    }
    return {{(*range)[0], (*range)[1], (*range)[2]}};
  }
  return {};
}
template <class Vector3>
optional<Vector3>
configuration::get_optional_vector3(const string_path& key,
                                    const translation_map<double>& tr) const {
  auto range = get_optional_vector(key, tr);
  if (range) {
    if (range->size() != 3) {
      throw translation_error(key, stringify(*range));
    }
    return {{(*range)[0], (*range)[1], (*range)[2]}};
  }
  return {};
}
template <class Vector3>
Vector3 configuration::get_vector3(const string_path& key) const {
  auto range = get_optional_vector3<Vector3>(key);
  if (!range) {
    throw key_error(key);
  }
  return *range;
}
template <class Vector3>
Vector3 configuration::get_vector3(const string_path& key,
                                   const translation_map<double>& tr) const {
  auto range = get_optional_vector3<Vector3>(key, tr);
  if (!range) {
    throw key_error(key);
  }
  return *range;
}
// and 4-Vector version
template <class Vector4>
optional<Vector4>
configuration::get_optional_vector4(const string_path& key) const {
  auto range = get_optional_vector<double>(key);
  if (range) {
    if (range->size() != 4) {
      throw translation_error(key, stringify(*range));
    }
    return {{(*range)[0], (*range)[1], (*range)[2], (*range)[3]}};
  }
  return {};
}
template <class Vector4>
optional<Vector4>
configuration::get_optional_vector4(const string_path& key,
                                    const translation_map<double>& tr) const {
  auto range = get_optional_vector(key, tr);
  if (range) {
    if (range->size() != 4) {
      throw translation_error(key, stringify(*range));
    }
    return {{(*range)[0], (*range)[1], (*range)[2], (*range)[3]}};
  }
  return {};
}
template <class Vector4>
Vector4 configuration::get_vector4(const string_path& key) const {
  auto range = get_optional_vector4<Vector4>(key);
  if (!range) {
    throw key_error(key);
  }
  return *range;
}
template <class Vector4>
Vector4 configuration::get_vector4(const string_path& key,
                                   const translation_map<double>& tr) const {
  auto range = get_optional_vector4<Vector4>(key, tr);
  if (!range) {
    throw key_error(key);
  }
  return *range;
}
// and bitpattern versiosn
template <class T>
optional<T>
configuration::get_optional_bitpattern(const string_path& key) const {
  optional<std::vector<T>> vec{get_optional_vector<T>(key)};
  optional<T> pattern;
  if (vec) {
    pattern.reset(static_cast<T>(0));
    for (const auto& val : *vec) {
      *pattern = static_cast<T>(*pattern | val);
    }
  }
  return pattern;
}
template <class T>
optional<T>
configuration::get_optional_bitpattern(const string_path& key,
                                       const translation_map<T>& tr) const {
  optional<std::vector<T>> vec{get_optional_vector(key, tr)};
  optional<T> pattern;
  if (vec) {
    pattern.reset(static_cast<T>(0));
    for (const auto& val : *vec) {
      *pattern = static_cast<T>(*pattern | val);
    }
  }
  return pattern;
}
template <class T>
T configuration::get_bitpattern(const string_path& key) const {
  auto s = get_optional_bitpattern<T>(key);
  if (!s) {
    throw key_error(key);
  }
  return *s;
}
template <class T>
T configuration::get_bitpattern(const string_path& key,
                                const translation_map<T>& tr) const {
  auto s = get_optional_bitpattern<T>(key, tr);
  if (!s) {
    throw key_error(key);
  }
  return *s;
}
// and "range" (interval) version
template <class T>
optional<interval<T>>
configuration::get_optional_range(const string_path& key) const {
  auto range = get_optional_vector<T>(key);
  if (range) {
    if (range->size() != 2) {
      throw translation_error(key, stringify(*range));
    }
    return {{(*range)[0], (*range)[1]}};
  }
  return {};
}
template <class T>
optional<interval<T>>
configuration::get_optional_range(const string_path& key,
                                  const translation_map<T>& tr) const {
  auto range = get_optional_vector(key, tr);
  if (range) {
    if (range->size() != 2) {
      throw translation_error(key, stringify(*range));
    }
    return {{(*range)[0], (*range)[1]}};
  }
  return {};
}
template <class T>
interval<T> configuration::get_range(const string_path& key) const {
  auto range = get_optional_range<T>(key);
  if (!range) {
    throw key_error(key);
  }
  return *range;
}
template <class T>
interval<T> configuration::get_range(const string_path& key,
                                     const translation_map<T>& tr) const {
  auto range = get_optional_range(key, tr);
  if (!range) {
    throw key_error(key);
  }
  return *range;
}

// =============================================================================
// Implementation: configuration errors
// =============================================================================

// configuration_translation_error<T> impl
template <class T>
configuration_translation_error
configuration::translation_error(const string_path& key,
                                 const std::string& value,
                                 const translation_map<T>& tr) const {
  return {key, value, tr, settings_path_, defaults_path_};
}

// "manual" translation (private)
template <class T>
T configuration::translate(const string_path& key, const std::string& val,
                           const translation_map<T>& tr) const {
  try {
    return tr.at(val);
  } catch (std::out_of_range) {
    throw translation_error(key, val, tr);
  }
}

// further configuration_translation_error implementation
template <class T>
configuration_translation_error::configuration_translation_error(
    const string_path& key, const std::string& value,
    const translation_map<T>& tr, const string_path& settings_path,
    const string_path& defaults_path)
    : configuration_error(
          "Unable to translate value '" + value + "' for key '" + key.str() +
              "' (in '" + settings_path.str() + "' or '" + defaults_path.str() +
              "' -- allowed values: '" +
              stringify(tr, "', '",
                        [](const typename translation_map<T>::value_type& el) {
                          return el.first;
                        }) +
              "')",
          "configuration_translation_error") {}

} // ns pcsim

#endif
