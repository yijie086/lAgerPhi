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

#ifndef LAGER_FRAMEWORK_LOADED
#define LAGER_FRAMEWORK_LOADED

#include <functional>
#include <string>
#include <vector>

#include <lager/core/configuration.hh>
#include <lager/core/logo.hh>

#include <boost/program_options.hpp>

namespace lager {
// error prototypes
class framework_error;
class framework_help;
class framework_file_error;
class framework_parse_error;
} // ns lager

// =============================================================================
// analyser::framework
//
// If you want to use the framework, call MAKE_LAGER_FRAMEWORK(function) at
// the end of your main source file. "function" is the analysis function you
// want to call. The analysis function should take the following inputs:
//    * conf: configuration object made from the input configuration file
//    * output: the base name for the output files.
//  Note:
//    * the framework takes care of all fancy exception handling
// =============================================================================
#define MAKE_LAGER_FRAMEWORK(function)                                         \
  int main(int argc, char* argv[]) {                                           \
    try {                                                                      \
      lager::print_logo();                                                     \
      lager::framework lager{argc, argv, (function)};                          \
      lager.run();                                                             \
      return 0;                                                                \
    } catch (...) {                                                            \
      return 1;                                                                \
    }                                                                          \
  }

namespace lager {
class framework {
public:
  using lager_function_type =
      std::function<int(const configuration& conf, const std::string& output)>;

  // setup the analysis framework
  framework(int argc, char* argv[], lager_function_type lager_function);
  ~framework();
  // run the analyzis framework
  int run() const;

private:
  boost::program_options::variables_map parse_arguments(int argc,
                                                        char* argv[]) const;

  // get our settings ptree used to initialize the configuration
  ptree get_settings() const;

  // suppress the ROOT signal handler
  int root_suppress_signals() const;

  // get an option from the command line, use the configuration file as fallback
  template <class T> T get_option(const std::string& key);

  int dummy_; // dummy used to ensure root_suppress_signals is run before
              // everything else
  boost::program_options::variables_map args_;
  configuration conf_;
  std::string output_;
  std::ofstream log_file_;
  lager_function_type lager_function_;
};
} // ns lager

// =============================================================================
// Definition: exceptions
// =============================================================================
namespace lager {
class framework_error : public lager::exception {
public:
  framework_error(const std::string& msg,
                  const std::string& type = "framework_error")
      : lager::exception{msg, type} {}
};
class framework_help : public framework_error {
public:
  framework_help(const std::string& program,
                 const boost::program_options::options_description& opts);

private:
  std::string
  message(const std::string& program,
          const boost::program_options::options_description& opts) const;
};
class framework_file_error : public framework_error {
public:
  framework_file_error(const std::string& file);
};
class framework_parse_error : public framework_error {
public:
  framework_parse_error(const std::string& file);
};
} // ns lager

// =============================================================================
// Implementation: framework
// =============================================================================
namespace lager {
// get an option from the command line, use the configuration file as fallback
template <class T> T framework::get_option(const std::string& key) {
  T val = 0;
  auto val_opt = conf_.get_optional<T>(key);
  if (args_.count(key)) {
    val = args_[key].as<T>();
    conf_.set(key, val);
  } else if (val_opt) {
    val = *val_opt;
  } else {
    throw framework_error(
        "Ensure that '" + key +
        "' is set on the command line or in the configuration file");
  }
  return val;
}
} // ns lager

#endif
