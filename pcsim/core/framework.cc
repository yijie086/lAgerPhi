#include "framework.hh"

#include <exception>
#include <cstdlib>

#include <pcsim/core/exception.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/core/stringify.hh>

#include <TSystem.h>

#include <boost/exception/diagnostic_information.hpp>
#include <boost/filesystem.hpp>

// useful aliases
namespace fs = boost::filesystem;
namespace po = boost::program_options;
// utility functions (unnamed namespace)
namespace {
void file_exists(const std::string& file) {
  if (!fs::exists(file)) {
    throw pcsim::framework_file_error{file};
  }
}
} // ns unnamed

namespace pcsim {
// =============================================================================
// framework constructor: suppress ROOT signals, parse command line arguments
// and provide for error handling
// =============================================================================
framework::framework(int argc, char* argv[],
                     pcsim_function_type pcsim_function)
    : pcsim_function_{pcsim_function} {
  try {
    LOG_INFO("pcsim", "Starting pcsim framework");
    // suppress ROOT signal handler
    root_suppress_signals();
    // parse the command line
    args_ = parse_arguments(argc, argv);
    // configuration
    settings_ = get_settings();
    LOG_INFO("pcsim", "Configuration file: " + args_["conf"].as<std::string>());
    // run and events
    int run = args_["run"].as<int>();
    int events = args_["events"].as<int>();
    settings_.put("mc.run", run);
    settings_.put("mc.events", events);
    // output root
    output_ = args_["out"].as<std::string>();
    char suffix[1024];
    sprintf(suffix, ".run%06i-%i", run, events);
    output_ += suffix;
    LOG_INFO("pcsim", "Output files will be written to: " + output_ + ".*");
  } catch (const framework_help& h) {
    std::cerr << h.what() << std::endl;
    exit(0);
  } catch (const pcsim::exception& e) {
    LOG_ERROR(e.type(), e.what());
    LOG_ERROR(e.type(), "Run with -h for help.");
    throw e;
  } catch (const boost::exception& e) {
    LOG_CRITICAL("boost::exception", boost::diagnostic_information(e));
    LOG_CRITICAL("boost::exception", "Unhandled boost exception");
    LOG_CRITICAL("boost::exception", "Please contact developer for support.");
    throw pcsim::exception("Unhandled boost exception", "boost::exception");
  } catch (const std::exception& e) {
    LOG_CRITICAL("std::exception", e.what());
    LOG_CRITICAL("std::exception", "Unhandled standard exception");
    LOG_CRITICAL("std::exception", "Please contact developer for support.");
    throw pcsim::exception("Unhandled standard exception", "std::exception");
  }
}
int framework::run() const {
  try{
    LOG_INFO("pcsim", "Starting event generator...");
    int ret = pcsim_function_(settings_, output_);
    LOG_INFO("pcsim", "Finished.");
    return ret;
  } catch (const pcsim::exception& e) {
    LOG_ERROR(e.type(), e.what());
    LOG_ERROR(e.type(), "Run with -h for help.");
    throw e;
  } catch (const boost::exception& e) {
    LOG_CRITICAL("boost::exception", boost::diagnostic_information(e));
    LOG_CRITICAL("boost::exception", "Unhandled boost exception");
    LOG_CRITICAL("boost::exception", "Please contact developer for support.");
    throw pcsim::exception("Unhandled boost exception", "boost::exception");
  } catch (const std::exception& e) {
    LOG_CRITICAL("std::exception", e.what());
    LOG_CRITICAL("std::exception", "Unhandled standard exception");
    LOG_CRITICAL("std::exception", "Please contact developer for support.");
    throw pcsim::exception("Unhandled standard exception", "std::exception");
  }
}
} // ns pcsim

// =============================================================================
// framework private utility functions 
// =============================================================================
namespace pcsim {
// =============================================================================
// Implementation: framework::parse_arguments
// Also sets the verbosity level to what was requested
// =============================================================================
po::variables_map framework::parse_arguments(int argc, char* argv[]) const {
  po::variables_map args;
  try {
    po::options_description opts_visible{"Allowed options"};
    opts_visible.add_options()("help,h", "Produce help message")(
        "conf,c", po::value<std::string>()->required()->notifier(file_exists),
        "Configuration JSON file")("run,r", po::value<int>()->required(),
                                   "Run number (also the random seed)")(
        "events,e", po::value<int>()->required(),
        "Number of events to generate")(
        "verb,v", po::value<unsigned>()->default_value(
                      static_cast<unsigned>(log_level::INFO)),
        "Verbosity level (0 -> 7; 0: silent, 4: default, 5: debug)")(
        "out,o", po::value<std::string>()->required(), "Output file name root");
    po::options_description opts_flags;
    opts_flags.add(opts_visible);
    po::positional_options_description opts_positional;

    po::store(po::command_line_parser(argc, argv)
                  .options(opts_flags)
                  .positional(opts_positional)
                  .run(),
              args);

    // help message requested? (BEFORE notify!)
    if (args.count("help")) {
      throw framework_help{argv[0], opts_visible};
    }
    // do our actual processing
    po::notify(args);
    // set the verbosity level if requested
    if (args.count("verb")) {
      unsigned v{args["verb"].as<unsigned>()};
      LOG_INFO("pcsim", "Verbosity level: " + std::to_string(v));
      global::logger.set_level(v);
    }
    return args;
  } catch (const po::error& e) {
    throw framework_error{e.what()};
  }
  return args;
}
// =============================================================================
// Implementation: framework::get_settings
// =============================================================================
ptree framework::get_settings() const {
  ptree settings;
  try {
    read_json(args_["conf"].as<std::string>(), settings);
  } catch (const boost::property_tree::ptree_error& e) {
    LOG_ERROR("framework_parse_error", e.what());
    throw framework_parse_error{args_["conf"].as<std::string>()};
  }
  return settings;
}
// =============================================================================
// Implementation: framework::root_suppress_signals
// Suppress the ROOT signal handlers, as they can cause undefined behavior and
// interfere with debugging
// =============================================================================
void framework::root_suppress_signals() const {
  gSystem->ResetSignal(kSigChild);
  gSystem->ResetSignal(kSigBus);
  gSystem->ResetSignal(kSigSegmentationViolation);
  gSystem->ResetSignal(kSigIllegalInstruction);
  gSystem->ResetSignal(kSigSystem);
  gSystem->ResetSignal(kSigPipe);
  gSystem->ResetSignal(kSigAlarm);
  gSystem->ResetSignal(kSigUrgent);
  gSystem->ResetSignal(kSigFloatingException);
  gSystem->ResetSignal(kSigWindowChanged);
}
} // ns pcsim

// =============================================================================
// Implementation: Exceptions
// =============================================================================
namespace pcsim {
framework_file_error::framework_file_error(const std::string& file)
    : framework_error{"No such file or directory: " + file,
                      "framework_file_error"} {}
framework_parse_error::framework_parse_error(const std::string& file)
    : framework_error{"Failed to parse: " + file, "framework_parse_error"} {}
framework_help::framework_help(const std::string& program,
                               const po::options_description& opts)
    : framework_error{message(program, opts), "help"} {}
std::string framework_help::message(const std::string& program,
                                    const po::options_description& opts) const {
  std::stringstream ss;
  ss << "\nUsage: " << program << " [options]\n" << opts
     << "\n";
  return ss.str();
}
} // ns pcsim
