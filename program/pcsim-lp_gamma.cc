#include <TFile.h>
#include <TRandom3.h>
#include <memory>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/framework.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/core/progress_meter.hh>
#include <pcsim/gen/lp_gamma_event.hh>
#include <pcsim/gen/lp_gamma_generator.hh>

using namespace pcsim;

// util function
// TODO add to core
std::string to_string_exp(double d) {
  std::stringstream ss;
  ss << std::scientific << d;
  return ss.str();
}

int run_mc(const configuration& cf, const std::string& output) {
  LOG_INFO("pcsim-lp_gamma", "Initializing PCSIM for lp-gamma processes");

  // get RNG
  LOG_INFO("pcsim-lp_gamma",
            "Initializing the RNG with seed " + cf.get<std::string>("run"));
  std::shared_ptr<TRandom> r {std::make_shared<TRandom3>()};
  r->SetSeed(cf.get<int>("run"));

  // make output file and buffer
  LOG_INFO("pcsim-lp_gamma", "Initializing the output buffer");
  std::shared_ptr<TFile> ofile{
      std::make_shared<TFile>((output + ".root").c_str(), "recreate")};
  lp_gamma_out evbuf{ofile, "lp_gamma_event"};

  // get event generator
  LOG_INFO("pcsim-lp_gamma", "Initializing the event generator");
  lp_gamma_generator gen{cf, "generator", r};

  // number of requested events:
  const size_t events = cf.get<size_t>("events");
  progress_meter progress{events};

  // loop over events
  LOG_INFO("pcsim-lp_gamma", "Starting the main generation loop");
  while (gen.n_events() < events) {
    evbuf.push(gen.generate());
    progress.update();
  }
  LOG_INFO("pcsim", "Event generation complete");
  LOG_INFO("pcsim", "Total number of generated events: " +
                        std::to_string(gen.n_events()));
  LOG_INFO("pcsim",
           "Total cross section [nb]: " + to_string_exp(gen.cross_section()));

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
