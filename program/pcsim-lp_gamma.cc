#include <TFile.h>
#include <TRandom3.h>
#include <memory>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/framework.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/core/progress_meter.hh>
#include <pcsim/gen/lp_gamma_event.hh>
#include <pcsim/gen/lp_gamma_generator.hh>

// TODO fix this
#include <pcsim/gen/beam/photon_gen.hh>
#include <pcsim/gen/beam/primary_gen.hh>
#include <pcsim/gen/lp_gamma/brodsky_2vmX.hh>
#include <pcsim/gen/lp_gamma/gaussian_1X.hh>
#include <pcsim/proc/detector/solid.hh>
#include <pcsim/proc/detector/jleic.hh>
#include <pcsim/proc/detector/null.hh>
// TODO

using namespace pcsim;

// util function
// TODO add to core
std::string to_string_exp(double d) {
  std::stringstream ss;
  ss << std::scientific << d;
  return ss.str();
}

int run_mc(const configuration& cf, const std::string& output) {

  // TODO fix this
  FACTORY_REGISTER2(lp_gamma::generator, lp_gamma::brodsky_2vmX,
                    "brodsky_2vmX");
  FACTORY_REGISTER2(lp_gamma::generator, lp_gamma::gaussian_qpq,
                    "gaussian_1qpq");
  FACTORY_REGISTER2(beam::primary_generator, beam::beam, "primary");

  FACTORY_REGISTER2(beam::photon_generator, beam::bremsstrahlung,
                    "bremsstrahlung");
  FACTORY_REGISTER2(beam::photon_generator, beam::vphoton, "vphoton");
  FACTORY_REGISTER2(detector::detector, detector::solid, "solid");
  FACTORY_REGISTER2(detector::detector, detector::jleic, "jleic");
  FACTORY_REGISTER2(detector::detector, detector::null, "4pi");
  // TODO

  LOG_INFO("pcsim-lp_gamma", "Initializing PCSIM for lp-gamma processes");

  // get RNG
  LOG_INFO("pcsim-lp_gamma",
           "Initializing the RNG with seed " + cf.get<std::string>("run"));
  std::shared_ptr<TRandom> r{std::make_shared<TRandom3>()};
  r->SetSeed(cf.get<int>("run"));

  // make output file and buffer
  LOG_INFO("pcsim-lp_gamma", "Initializing the output buffer");
  std::shared_ptr<TFile> ofile{
      std::make_shared<TFile>((output + ".root").c_str(), "recreate")};
  lp_gamma_out evbuf{ofile, "lp_gamma_event"};

  // get event generator
  LOG_INFO("pcsim-lp_gamma", "Initializing the event generator");
  lp_gamma_generator gen{cf, "generator", r};

  // init the progress meter with number of requested events
  progress_meter progress{static_cast<size_t>(gen.n_requested())};

  // loop over events
  LOG_INFO("pcsim-lp_gamma", "Starting the main generation loop");
  while (!gen.finished()) {
    evbuf.push(gen.generate());
    progress.update(gen.n_events(), gen.n_requested());
  }
  LOG_INFO("pcsim", "Event generation complete");
  LOG_INFO("pcsim", "Total number of generated events: " +
                        std::to_string(gen.n_events()));
  LOG_INFO("pcsim",
           "Total cross section [nb]: " + to_string_exp(gen.cross_section()));

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
