#include "configuration.hh"

#include <pcsim/core/logger.hh>

using boost::property_tree::ptree_bad_path;
using boost::property_tree::ptree_error;

// =============================================================================
// class configuration
// =============================================================================
namespace pcsim {
configuration::configuration(const ptree& settings, const string_path& path)
    : settings_path_{path} {
  load(settings);
}
// load
void configuration::load(const ptree& in_conf) {
  try {

    LOG_DEBUG(settings_path_.str(), "Loading settings");
    auto set = in_conf.get_child_optional(settings_path_);
    if (!set) {
      throw path_error(settings_path_.str());
    }
    settings_ = *set;

    auto type = settings_.get_optional<std::string>(TYPE_KEY);
    if (!type) {
      throw configuration_error("type descriptor '" + std::string(TYPE_KEY) +
                                "' has to be set in " + settings_path_.str());
    }
    defaults_path_ /= DEFAULTS;
    defaults_path_ /= *type;
    auto def = in_conf.get_child_optional(defaults_path_);
    if (def) {
      LOG_JUNK(settings_path_.str(), *type + " defaults found.");
      defaults_ = *def;
    } else {
      LOG_JUNK(settings_path_.str(),
               "No default settings provided for this module.");
    }
  } catch (ptree_error& e) {
    // shouldn't happen
    throw configuration_error("Processing error", e.what());
  } 
}
// save
void configuration::save(ptree& out_conf) const {
  out_conf.put_child(settings_path_, settings_);
  LOG_DEBUG(settings_path_.str(), "Settings saved.");
  if (!out_conf.get_child_optional(defaults_path_)) {
    out_conf.put_child(defaults_path_, defaults_);
    LOG_INFO(defaults_path_.str(), "Settings saved.");
  }
}

configuration_path_error
configuration::path_error(const string_path& path) const {
  return {path};
}
configuration_key_error configuration::key_error(const string_path& key) const {
  return {key, settings_path_, defaults_path_};
}
configuration_value_error
configuration::value_error(const string_path& key,
                           const std::string& value) const {
  return {key, value, settings_path_, defaults_path_};
}
configuration_translation_error
configuration::translation_error(const string_path& key,
                                 const std::string& value) const {
  return {key, value, settings_path_, defaults_path_};
}
} // ns pcsim
// =============================================================================
// exceptions
// =============================================================================
namespace pcsim {
// path error
configuration_path_error::configuration_path_error(const string_path& path)
    : configuration_error{"Cannot find the configuration path '" +
                              path.str() + "'",
                          "configuration_path_error"} {}
// key error
configuration_key_error::configuration_key_error(
    const string_path& key, const string_path& settings_path,
    const string_path& defaults_path)
    : configuration_error{"Cannot find '" + key.str() + "' (in '" +
                              settings_path.str() + "' or '" +
                              defaults_path.str() + "')",
                          "configuration_key_error"} {}
// value error
configuration_value_error::configuration_value_error(
    const string_path& key, const std::string& value,
    const string_path& settings_path, const string_path& defaults_path)
    : configuration_error{"Invalid value '" + value + "' for key '" +
                              key.str() + "' (in '" + settings_path.str() +
                              "' or '" + defaults_path.str() + "')",
                          "configuration_value_error"} {}

// translation error
configuration_translation_error::configuration_translation_error(
    const string_path& key, const std::string& value,
    const string_path& settings_path, const string_path& defaults_path)
    : configuration_error{"Unable to translate value '" + value +
                              "' for key '" + key.str() + "' (in '" +
                              settings_path.str() + "' or '" +
                              defaults_path.str() + "')",
                          "configuration_translation_error"} {}
} // ns pcsim

// =============================================================================
// utility functions to format object names and titles
// =============================================================================
namespace pcsim {
std::string format_name(std::string context, const std::string& name) {
  if (name.size() && context.size()) {
    context += NAME_SEPARATOR;
    context += name;
    return context;
  } else if (context.size()) {
    return context;
  }
  return name;
}
std::string format_title(const std::string& context, std::string title) {
  if (title.size() && context.size()) {
    title += TITLE_SEPARATOR;
    title += context;
  } else if (context.size()) {
    return context;
  }
  return title;
}
} // ns pcsim
