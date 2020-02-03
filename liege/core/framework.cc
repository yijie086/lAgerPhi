#include "framework.hh"
#include "version.hh"

#include <cstdlib>
#include <exception>

#include <liege/core/configuration.hh>
#include <liege/core/exception.hh>
#include <liege/core/logger.hh>
#include <liege/core/stringify.hh>

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
    throw liege::framework_file_error{file};
  }
}
} // namespace

namespace liege {
// =============================================================================
// framework constructor: suppress ROOT signals, parse command line arguments
// and provide for error handling
// =============================================================================
framework::framework(int argc, char* argv[],
                     liege_function_type liege_function) try
    // suppress ROOT signal handler
    : dummy_ {
  root_suppress_signals()
}
// parse the command line
, args_{parse_arguments(argc, argv)} // configuration
,
    conf_{get_settings(), "mc"} // framework function
,
    liege_function_{liege_function} {
  // talk to the user
  LOG_INFO("liege", "Starting liege framework");
  LOG_INFO("liege", "Configuration file: " + args_["conf"].as<std::string>());

  // get our run info and luminosity/number of events to be generated
  const int run = get_option<int>("run");
  auto lumi = conf_.get_optional<double>("lumi");
  int events = -1;
  if (!lumi) {
    events = get_option<int>("events");
  }

  // output file name

  // path
  output_ = args_["out"].as<std::string>();
  if (output_.back() != '/') {
    output_ += '/';
  }

  // what MC program are we running?
  output_ += conf_.get<std::string>("type");

  // add generator name and acceptance simulation
  output_ += "." + conf_.get<std::string>("generator/type");

  auto detector = conf_.get<std::string>("detector/type");
  if (detector.size() > 0 && detector != "none") {
    output_ += "." + conf_.get<std::string>("detector/type");
  } else {
    output_ += ".4pi";
  }

  // add optional tag
  auto tag = conf_.get_optional<std::string>("tag");
  if (tag && tag->size()) {
    output_ += "." + *tag;
  }
  // add the run info and lumi or number of generated events
  char info[1024];
  if (lumi) {
    sprintf(info, ".run%05i-lumi%f", run, *lumi);
  } else {
    sprintf(info, ".run%05i-%i", run, events);
  }
  output_ += info;

  // communicate file name to user
  LOG_INFO("liege", "Output files will be written to: " + output_ + ".*");

  // write the a copy of the configuration file to the output
  ptree settings;
  conf_.save(settings);
  write_json(output_ + ".json", settings);

  // redirect logger to use the log file
  LOG_INFO("liege", "Redirecting logger to: " + output_ + ".log");
  log_file_.open(output_ + ".log");
  global::logger.set_output(log_file_);

  // Communicate LIEGE commit hash for this run
  LOG_INFO("liege", "LIEGE version: " + get_liege_version());
}
catch (const framework_help& h) {
  std::cerr << h.what() << std::endl;
  exit(0);
}
catch (const liege::exception& e) {
  LOG_ERROR(e.type(), e.what());
  LOG_ERROR(e.type(), "Run with -h for help.");
  throw e;
}
catch (const boost::exception& e) {
  LOG_CRITICAL("boost::exception", boost::diagnostic_information(e));
  LOG_CRITICAL("boost::exception", "Unhandled boost exception");
  LOG_CRITICAL("boost::exception", "Please contact developer for support.");
  throw liege::exception("Unhandled boost exception", "boost::exception");
}
catch (const std::exception& e) {
  LOG_CRITICAL("std::exception", e.what());
  LOG_CRITICAL("std::exception", "Unhandled standard exception");
  LOG_CRITICAL("std::exception", "Please contact developer for support.");
  throw liege::exception("Unhandled standard exception", "std::exception");
}

framework::~framework() {
  // redirect output stream back to the standard output
  global::logger.set_output(std::cout);
}

int framework::run() const {
  try {
    LOG_INFO("liege", "Starting event generator...");
    int ret = liege_function_(conf_, output_);
    LOG_INFO("liege", "Finished.");
    return ret;
  } catch (const liege::exception& e) {
    LOG_ERROR(e.type(), e.what());
    LOG_ERROR(e.type(), "Run with -h for help.");
    throw e;
  } catch (const boost::exception& e) {
    LOG_CRITICAL("boost::exception", boost::diagnostic_information(e));
    LOG_CRITICAL("boost::exception", "Unhandled boost exception");
    LOG_CRITICAL("boost::exception", "Please contact developer for support.");
    throw liege::exception("Unhandled boost exception", "boost::exception");
  } catch (const std::exception& e) {
    LOG_CRITICAL("std::exception", e.what());
    LOG_CRITICAL("std::exception", "Unhandled standard exception");
    LOG_CRITICAL("std::exception", "Please contact developer for support.");
    throw liege::exception("Unhandled standard exception", "std::exception");
  }
}

} // namespace liege

// =============================================================================
// framework private utility functions
// =============================================================================
namespace liege {
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
        "Configuration JSON file")("run,r", po::value<int>(),
                                   "Run number (also the random seed)")(
        "events,e", po::value<int>(), "Number of events to generate")(
        "verb,v",
        po::value<unsigned>()->default_value(
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
      LOG_INFO("liege", "Verbosity level: " + std::to_string(v));
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
int framework::root_suppress_signals() const {
  return 0; // this somehow does not work anymore after the newest ROOT update
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
  return 0;
}
} // namespace liege

// =============================================================================
// Implementation: Exceptions
// =============================================================================
namespace liege {
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
  ss << "\nUsage: " << program << " [options]\n" << opts << "\n";
  return ss.str();
}
} // namespace liege
