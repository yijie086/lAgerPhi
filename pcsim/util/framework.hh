#ifndef PYTHIA6M_FRAMEWORK_LOADED
#define PYTHIA6M_FRAMEWORK_LOADED

#include <functional>
#include <vector>
#include <string>

#include <pythia6m/util/configuration.hh>

#include <boost/program_options.hpp>

namespace pythia6m {
// error prototypes
class framework_error;
class framework_help;
class framework_file_error;
class framework_parse_error;
} // ns pythia6m

// =============================================================================
// analyser::framework
//
// If you want to use the framework, call MAKE_PYTHIA6M_FRAMEWORK(function) at
// the end of your main source file. "function" is the analysis function you
// want to call. The analysis function should take the following inputs:
//    * settings: ptree made from the input configuration file
//    * output: the base name for the output files.
//  Note:
//    * the framework takes care of all fancy exception handling
// =============================================================================
#define MAKE_PYTHIA6M_FRAMEWORK(function)                                      \
  int main(int argc, char* argv[]) {                                           \
    try {                                                                      \
      pythia6m::framework pythia6m{argc, argv, (function)};                    \
      pythia6m.run();                                                          \
      return 0;                                                                \
    } catch (...) {                                                            \
      return 1;                                                                \
    }                                                                          \
  }

namespace pythia6m {
class framework {
public:
  using pythia6m_function_type =
      std::function<int(const ptree& settings, const std::string& output)>;

  // setup the analysis framework
  framework(int argc, char* argv[], pythia6m_function_type pythia6m_function);
  // run the analyzis framework
  int run() const;

private:
  boost::program_options::variables_map parse_arguments(int argc,
                                                        char* argv[]) const;
  ptree get_settings() const;
  void root_suppress_signals() const;

  pythia6m_function_type pythia6m_function_;
  boost::program_options::variables_map args_;
  std::string output_;
  ptree settings_;
};
} // ns pythia6m

// =============================================================================
// Definition: exceptions
// =============================================================================
namespace pythia6m {
class framework_error : public pythia6m::exception {
public:
  framework_error(const std::string& msg,
                  const std::string& type = "framework_error")
      : pythia6m::exception{msg, type} {}
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
} // ns pythia6m

#endif
