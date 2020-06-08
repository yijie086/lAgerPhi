// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2020 Sylvester Joosten <sjoosten@anl.gov>
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

#include <lager/core/configuration.hh>
#include <lager/core/framework.hh>
#include <lager/core/logger.hh>
#include <lager/core/progress_meter.hh>
#include <lager/gen/lA_event.hh>
#include <lager/gen/lA_generator.hh>

#include <HepMC3/WriterAscii.h>
#include <TFile.h>
#include <TRandom3.h>
#include <fstream>
#include <memory>

// TODO fix this
#include <lager/gen/initial/beam_gen.hh>
#include <lager/gen/initial/photon_gen.hh>
#include <lager/gen/initial/target_gen.hh>
#include <lager/gen/initial/vertex_gen.hh>
#include <lager/gen/lA/brodsky_2vmX.hh>
#include <lager/gen/lA/lee_4He_jpsi_grid.hh>
#include <lager/gen/lA/oleksii_2vmp.hh>
#include <lager/gen/lA/oleksii_jpsi_bh.hh>
#include <lager/gen/lA/resonance_qpq.hh>
#include <lager/proc/detector/composite.hh>
#include <lager/proc/detector/cone.hh>
#include <lager/proc/detector/null.hh>
#include <lager/proc/detector/spectrometer.hh>
// TODO

using namespace lager;

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
  FACTORY_REGISTER2(lA::generator, lA::lee_4He_jpsi_grid, "lee_4He_jpsi_grid");
  FACTORY_REGISTER2(lA::generator, lA::brodsky_2vmX, "brodsky_2vmX");
  FACTORY_REGISTER2(lA::generator, lA::oleksii_2vmp, "oleksii_2vmp");
  FACTORY_REGISTER2(lA::generator, lA::oleksii_jpsi_bh, "oleksii_jpsi_bh");
  FACTORY_REGISTER2(lA::generator, lA::resonance_qpq, "resonance_1qpq");
  FACTORY_REGISTER2(initial::vertex_generator, initial::origin_vertex,
                    "origin");
  FACTORY_REGISTER2(initial::vertex_generator, initial::linear_vertex,
                    "linear");
  FACTORY_REGISTER2(initial::beam_generator, initial::constant_beam,
                    "constant");
  FACTORY_REGISTER2(initial::target_generator, initial::primary_target,
                    "primary");
  FACTORY_REGISTER2(initial::target_generator, initial::fermi87, "fermi87");

  FACTORY_REGISTER2(initial::photon_generator, initial::no_photon, "no-photon");
  FACTORY_REGISTER2(initial::photon_generator, initial::bremsstrahlung,
                    "bremsstrahlung");
  FACTORY_REGISTER2(initial::photon_generator,
                    initial::bremsstrahlung_realistic_target,
                    "bremsstrahlung_realistic_target");
  FACTORY_REGISTER2(initial::photon_generator, initial::vphoton, "vphoton");
  FACTORY_REGISTER2(detector::detector, detector::null, "4pi");
  FACTORY_REGISTER2(detector::detector, detector::spectrometer, "spectrometer");
  FACTORY_REGISTER2(detector::detector, detector::cone, "cone");
  FACTORY_REGISTER2(detector::detector, detector::composite, "composite");
  // TODO

  LOG_INFO("lager", "Initializing LAGER for lp-gamma processes");

  // get RNG
  LOG_INFO("lager",
           "Initializing the RNG with seed " + cf.get<std::string>("run"));
  std::shared_ptr<TRandom> r{std::make_shared<TRandom3>()};
  r->SetSeed(cf.get<int>("run"));

  // make output file and buffer
  LOG_INFO("lager", "Initializing the output buffer");
  std::shared_ptr<TFile> ofile{
      std::make_shared<TFile>((output + ".root").c_str(), "recreate")};

  // check if we want hepmc output as well
  std::unique_ptr<HepMC3::WriterAscii> ohepmc;
  auto do_hepmc = cf.get_optional<bool>("output_hepmc");
  if (do_hepmc && *do_hepmc) {
    LOG_INFO("lager", "Also outputting text output for HepMC");
    ohepmc =
        std::make_unique<HepMC3::WriterAscii>((output + ".hepmc").c_str());
  }
  // check if we want gemc output as well
  std::unique_ptr<std::ofstream> ogemc;
  auto do_gemc = cf.get_optional<bool>("output_gemc");
  if (do_gemc && *do_gemc) {
    LOG_INFO("lager", "Also outputting text output for GEMC");
    ogemc = std::make_unique<std::ofstream>(output + ".gemc");
  }
  // check if we want simc, in similar vein
  std::unique_ptr<std::ofstream> osimc;
  auto do_simc = cf.get_optional<bool>("output_simc");
  if (do_simc && *do_simc) {
    LOG_INFO("lager", "Also outputting text output for SIMC");
    osimc = std::make_unique<std::ofstream>(output + ".simc");
  }

  lA_out evbuf{ofile, std::move(ohepmc), std::move(ogemc), std::move(osimc),
               "lAger"};
  // get event generator
  LOG_INFO("lager", "Initializing the event generator");
  lA_generator gen{cf, "generator", r};

  // init the progress meter with number of requested events
  progress_meter progress{static_cast<size_t>(gen.n_requested())};

  // loop over events
  LOG_INFO("lager", "Starting the main generation loop");
  while (!gen.finished()) {
    evbuf.push(gen.generate());
    progress.update(gen.n_events(), gen.n_requested());
  }

  LOG_INFO("lager", "Event generation complete");
  LOG_INFO("lager",
           "Total number of generated events: " +
               std::to_string(gen.n_events()));
  LOG_INFO("lager",
           "Total accepted cross section [nb]: " +
               to_string_exp(gen.cross_section()));
  LOG_INFO("lager",
           "Partial accepted cross section with BR [nb]: " +
               to_string_exp(gen.partial_cross_section()));
  LOG_INFO("lager",
           " --> Acceptance [%]: " + std::to_string(100 * gen.acceptance()));
  // write generation statistics to file as 1D histograms
  LOG_INFO("lager", "Writing generation statistics to output file");
  write_value_to_file(ofile, "weighted_cross_section",
                      gen.cross_section() * gen.n_events());
  write_value_to_file(ofile, "weighted_partial_cross_section",
                      gen.partial_cross_section() * gen.n_events());
  write_value_to_file(ofile, "n_events", gen.n_events());

  return 0;
}

MAKE_LAGER_FRAMEWORK(run_mc)
