#include <TFile.h>
#include <TRandom3.h>
#include <fstream>
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
#include <pcsim/gen/beam/vertex_gen.hh>
#include <pcsim/gen/lp_gamma/brodsky_2vmX.hh>
#include <pcsim/gen/lp_gamma/gaussian_1X.hh>
#include <pcsim/gen/lp_gamma/lee_4He_jpsi_grid.hh>
#include <pcsim/gen/lp_gamma/oleksii_2vmp.hh>
#include <pcsim/gen/lp_gamma/oleksii_jpsi_bh.hh>
#include <pcsim/proc/detector/barrel.hh>
#include <pcsim/proc/detector/composite.hh>
#include <pcsim/proc/detector/null.hh>
#include <pcsim/proc/detector/solid.hh>
#include <pcsim/proc/detector/spectrometer.hh>
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
  FACTORY_REGISTER2(lp_gamma::generator, lp_gamma::oleksii_2vmp,
                    "oleksii_2vmp");
  FACTORY_REGISTER2(lp_gamma::generator, lp_gamma::lee_4He_jpsi_grid,
                    "lee_4He_jpsi_grid");
  FACTORY_REGISTER2(lp_gamma::generator, lp_gamma::oleksii_jpsi_bh,
                    "oleksii_jpsi_bh");
  FACTORY_REGISTER2(lp_gamma::generator, lp_gamma::gaussian_qpq,
                    "gaussian_1qpq");
  FACTORY_REGISTER2(beam::vertex_generator, beam::origin_vertex, "origin");
  FACTORY_REGISTER2(beam::vertex_generator, beam::linear_vertex, "linear");
  FACTORY_REGISTER2(beam::primary_generator, beam::beam, "primary");

  FACTORY_REGISTER2(beam::photon_generator, beam::bremsstrahlung,
                    "bremsstrahlung");
  FACTORY_REGISTER2(beam::photon_generator,
                    beam::bremsstrahlung_realistic_target,
                    "bremsstrahlung_realistic_target");
  FACTORY_REGISTER2(beam::photon_generator, beam::vphoton, "vphoton");
  FACTORY_REGISTER2(detector::detector, detector::solid, "solid");
  FACTORY_REGISTER2(detector::detector, detector::null, "4pi");
  FACTORY_REGISTER2(detector::detector, detector::spectrometer, "spectrometer");
  FACTORY_REGISTER2(detector::detector, detector::barrel, "barrel");
  FACTORY_REGISTER2(detector::detector, detector::composite, "composite");
  // TODO

  LOG_INFO("pcsim", "Initializing PCSIM for lp-gamma processes");

  // get RNG
  LOG_INFO("pcsim",
           "Initializing the RNG with seed " + cf.get<std::string>("run"));
  std::shared_ptr<TRandom> r{std::make_shared<TRandom3>()};
  r->SetSeed(cf.get<int>("run"));

  // make output file and buffer
  LOG_INFO("pcsim", "Initializing the output buffer");
  std::shared_ptr<TFile> ofile{
      std::make_shared<TFile>((output + ".root").c_str(), "recreate")};

  // check if we want gemc output as well
  std::unique_ptr<std::ofstream> olund;
  auto do_gemc = cf.get_optional<bool>("output_gemc");
  if (do_gemc && *do_gemc) {
    LOG_INFO("pcsim", "Also outputting text output for GEMC");
    olund = std::make_unique<std::ofstream>(output + ".gemc.dat");
  }
  // check if we want simc, in similar vein
  std::unique_ptr<std::ofstream> osimc;
  auto do_simc = cf.get_optional<bool>("output_simc");
  if (do_simc && *do_simc) {
    LOG_INFO("pcsim", "Also outputting text output for SIMC");
    osimc = std::make_unique<std::ofstream>(output + ".simc.dat");
  }

  lp_gamma_out evbuf{ofile, std::move(olund), std::move(osimc),
                     "lp_gamma_event"};
  // get event generator
  LOG_INFO("pcsim", "Initializing the event generator");
  lp_gamma_generator gen{cf, "generator", r};

  // init the progress meter with number of requested events
  progress_meter progress{static_cast<size_t>(gen.n_requested())};

  // loop over events
  LOG_INFO("pcsim", "Starting the main generation loop");
  while (!gen.finished()) {
    evbuf.push(gen.generate());
    progress.update(gen.n_events(), gen.n_requested());
  }
  LOG_INFO("pcsim", "Event generation complete");
  LOG_INFO("pcsim", "Total number of generated events: " +
                        std::to_string(gen.n_events()));
  LOG_INFO("pcsim", "Total accepted cross section [nb]: " +
                        to_string_exp(gen.cross_section()));
  LOG_INFO("pcsim", "Partial accepted cross section with BR [nb]: " +
                        to_string_exp(gen.partial_cross_section()));
  LOG_INFO("pcsim",
           " --> Acceptance [%]: " + std::to_string(100 * gen.acceptance()));

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
