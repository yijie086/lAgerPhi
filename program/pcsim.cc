#include <pcsim/core/configuration.hh>
#include <pcsim/core/exception.hh>
#include <pcsim/core/framework.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/core/assert.hh>
#include <pcsim/gen/photon_gen.hh>
#include <memory>
#include <TRandom3.h>

using namespace pcsim;

int run_mc(const ptree& settings, const std::string& output) { 
  LOG_INFO("main", "PCSIM"); 

  // configuration
  configuration conf{settings, "mc"};
  const size_t events = conf.get<size_t>("events");

  // output file
  std::shared_ptr<TFile> ofile{
      std::make_shared<TFile>((output + ".root").c_str(), "recreate")};

  // init the RNG, use the run number as random seed
  std::shared_ptr<TRandom> rng {std::make_shared<TRandom3>()};
  rng->SetSeed(conf.get<int>("run"));

  // init the photon beam, add some diagnostic histos
  photon_gen pg(settings, "mc/photon_gen", rng);
  pg.add_histo(ofile, "Egamma", "Photon Beam Energy Spectrum",
               {"Egamma", [](const photon_beam& b) { return b.energy; }, 100,
                pg.range()});

  // generate our events
  for (size_t iev =0; iev < events; ++iev) {
    auto photon = pg.generate();
  }
}

MAKE_PCSIM_FRAMEWORK(run_mc)
