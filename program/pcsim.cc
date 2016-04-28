#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/exception.hh>
#include <pcsim/core/framework.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/core/progress_meter.hh>
#include <pcsim/gen/bremsstrahlung.hh>
#include <pcsim/gen/jpsi.hh>

using namespace pcsim;

// initialize the tree for J/Psi event data
void init_tree(TTree& t, gen::jpsi_event& evbuf) {
  t.Branch("xsec", &evbuf.xsec);
  t.Branch("flux", &evbuf.flux);
  t.Branch("weight", &evbuf.weight);
  t.Branch("s", &evbuf.s);
  t.Branch("t", &evbuf.t);
  t.Branch("tmin", &evbuf.tmin);
  t.Branch("tmax", &evbuf.tmax);
  t.Branch("W", &evbuf.W);
  t.Branch("beam", &evbuf.beam);
  t.Branch("target", &evbuf.target);
  t.Branch("recoil", &evbuf.recoil);
  t.Branch("jpsi", &evbuf.jpsi);
  t.Branch("positron", &evbuf.positron);
  t.Branch("electron", &evbuf.electron);
}

int run_mc(const ptree& settings, const std::string& output) { 
  LOG_INFO("main", "PCSIM"); 

  // configuration
  configuration conf{settings, "mc"};
  const size_t events = conf.get<size_t>("events");

  // output file
  std::shared_ptr<TFile> ofile{
      std::make_shared<TFile>((output + ".root").c_str(), "recreate")};
  ofile->cd();

  // output event tree and data buffer
  TTree* t = new TTree{"jpsi_event", "J/Psi Event Data"};
  gen::jpsi_event evbuf;
  init_tree(*t, evbuf);

  // init the RNG, use the run number as random seed
  std::shared_ptr<TRandom> rng {std::make_shared<TRandom3>()};
  rng->SetSeed(conf.get<int>("run"));

  // init the photon beam, add some diagnostic histos
  gen::bremsstrahlung photon_gen(settings, "mc/photon_gen", rng);
  photon_gen.add_histo(ofile, "Egamma", "Photon Energy",
                       {"Egamma",
                        [](const gen::photon_beam& b) { return b.energy; }, 100,
                        photon_gen.range()});

  // init the J/Psi generator, add some diagnostic histos
  gen::jpsi jpsi_gen(settings, "mc/jpsi_gen", rng);
  // progress meter
  progress_meter progress{events};

  // generate our events
  for (size_t iev = 0; iev < events; ++iev) {
    auto photon = photon_gen.generate();
    evbuf = jpsi_gen.generate(photon);
    if (!evbuf.good) {
      // rewind
      iev -= 1;
      continue;
    }
    t->Fill();
    progress.update();
  }

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
