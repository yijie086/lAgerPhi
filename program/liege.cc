#include <TFile.h>
#include <TRandom3.h>
#include <fstream>
#include <memory>
#include <liege/core/configuration.hh>
#include <liege/core/framework.hh>
#include <liege/core/logger.hh>
#include <liege/core/progress_meter.hh>
#include <liege/gen/lp_gamma_event.hh>
#include <liege/gen/lp_gamma_generator.hh>

// TODO fix this
#include <liege/gen/beam/photon_gen.hh>
#include <liege/gen/beam/primary_gen.hh>
#include <liege/gen/beam/vertex_gen.hh>
#include <liege/gen/lp_gamma/brodsky_2vmX.hh>
#include <liege/gen/lp_gamma/gaussian_1X.hh>
#include <liege/gen/lp_gamma/lee_4He_jpsi_grid.hh>
#include <liege/gen/lp_gamma/oleksii_2vmp.hh>
#include <liege/gen/lp_gamma/oleksii_jpsi_bh.hh>
#include <liege/proc/detector/barrel.hh>
#include <liege/proc/detector/composite.hh>
#include <liege/proc/detector/null.hh>
#include <liege/proc/detector/spectrometer.hh>
// TODO

using namespace liege;

// util function
// TODO add to core
std::string to_string_exp(double d) {
  std::stringstream ss;
  ss << std::scientific << d;
  return ss.str();
}

void write_value_to_file(std::shared_ptr<TFile> ofile, const std::string& name,
                         double value) {
  TH1D* tmp = new TH1D(name.c_str(), "", 1, 0, 1);
  tmp->SetBinContent(1, value);
  ofile->cd();
  tmp->Write();
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
  FACTORY_REGISTER2(detector::detector, detector::null, "4pi");
  FACTORY_REGISTER2(detector::detector, detector::spectrometer, "spectrometer");
  FACTORY_REGISTER2(detector::detector, detector::barrel, "barrel");
  FACTORY_REGISTER2(detector::detector, detector::composite, "composite");
  // TODO

  LOG_INFO("liege", "Initializing PCSIM for lp-gamma processes");

  // get RNG
  LOG_INFO("liege",
           "Initializing the RNG with seed " + cf.get<std::string>("run"));
  std::shared_ptr<TRandom> r{std::make_shared<TRandom3>()};
  r->SetSeed(cf.get<int>("run"));

  // make output file and buffer
  LOG_INFO("liege", "Initializing the output buffer");
  std::shared_ptr<TFile> ofile{
      std::make_shared<TFile>((output + ".root").c_str(), "recreate")};

  // check if we want gemc output as well
  std::unique_ptr<std::ofstream> olund;
  auto do_gemc = cf.get_optional<bool>("output_gemc");
  if (do_gemc && *do_gemc) {
    LOG_INFO("liege", "Also outputting text output for GEMC");
    olund = std::make_unique<std::ofstream>(output + ".gemc.dat");
  }
  // check if we want simc, in similar vein
  std::unique_ptr<std::ofstream> osimc;
  auto do_simc = cf.get_optional<bool>("output_simc");
  if (do_simc && *do_simc) {
    LOG_INFO("liege", "Also outputting text output for SIMC");
    osimc = std::make_unique<std::ofstream>(output + ".simc.dat");
  }

  lp_gamma_out evbuf{ofile, std::move(olund), std::move(osimc),
                     "lp_gamma_event"};
  // get event generator
  LOG_INFO("liege", "Initializing the event generator");
  lp_gamma_generator gen{cf, "generator", r};

  // init the progress meter with number of requested events
  progress_meter progress{static_cast<size_t>(gen.n_requested())};

  // loop over events
  LOG_INFO("liege", "Starting the main generation loop");
  while (!gen.finished()) {
    evbuf.push(gen.generate());
    progress.update(gen.n_events(), gen.n_requested());
  }

  LOG_INFO("liege", "Event generation complete");
  LOG_INFO("liege", "Total number of generated events: " +
                        std::to_string(gen.n_events()));
  LOG_INFO("liege", "Total accepted cross section [nb]: " +
                        to_string_exp(gen.cross_section()));
  LOG_INFO("liege", "Partial accepted cross section with BR [nb]: " +
                        to_string_exp(gen.partial_cross_section()));
  LOG_INFO("liege",
           " --> Acceptance [%]: " + std::to_string(100 * gen.acceptance()));
  // write generation statistics to file as 1D histograms
  LOG_INFO("liege", "Writing generation statistics to output file");
  write_value_to_file(ofile, "weighted_cross_section",
                      gen.cross_section() * gen.n_events());
  write_value_to_file(ofile, "weighted_partial_cross_section",
                      gen.partial_cross_section() * gen.n_events());
  write_value_to_file(ofile, "n_events", gen.n_events());

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
