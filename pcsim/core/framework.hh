#ifndef PCSIM_FRAMEWORK_LOADED
#define PCSIM_FRAMEWORK_LOADED

#include <functional>
#include <vector>
#include <string>

#include <pcsim/core/configuration.hh>

#include <boost/program_options.hpp>

namespace pcsim {
// error prototypes
class framework_error;
class framework_help;
class framework_file_error;
class framework_parse_error;
} // ns pcsim

// =============================================================================
// analyser::framework
//
// If you want to use the framework, call MAKE_PCSIM_FRAMEWORK(function) at
// the end of your main source file. "function" is the analysis function you
// want to call. The analysis function should take the following inputs:
//    * conf: configuration object made from the input configuration file
//    * output: the base name for the output files.
//  Note:
//    * the framework takes care of all fancy exception handling
// =============================================================================
#define MAKE_PCSIM_FRAMEWORK(function)                                         \
  int main(int argc, char* argv[]) {                                           \
    try {                                                                      \
      pcsim::framework pcsim{argc, argv, (function)};                          \
      pcsim.run();                                                             \
      return 0;                                                                \
    } catch (...) {                                                            \
      return 1;                                                                \
    }                                                                          \
  }

namespace pcsim {
class framework {
public:
  using pcsim_function_type =
      std::function<int(const configuration& conf, const std::string& output)>;

  // setup the analysis framework
  framework(int argc, char* argv[], pcsim_function_type pcsim_function);
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
  pcsim_function_type pcsim_function_;

};
} // ns pcsim

// =============================================================================
// Definition: exceptions
// =============================================================================
namespace pcsim {
class framework_error : public pcsim::exception {
public:
  framework_error(const std::string& msg,
                  const std::string& type = "framework_error")
      : pcsim::exception{msg, type} {}
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
} // ns pcsim


// =============================================================================
// Implementation: framework
// =============================================================================
namespace pcsim {
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
} // ns pcsim

#endif
