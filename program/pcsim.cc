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

// notes:
//  * For a give number of incoming electrons Ne, the number of good
//    bremsstrahlung photons Nb that for the process of interest is given by: 
//      Nb = Ne * flux * (evgen/ngamma)
//  * When binning using the cross section as event weight w_i, the differential
//    cross section in a bin is given by: 
//      dsigma = sum_i(w_i) / evgen / delta
//    with delta the bin width
struct mc_event {
  size_t event = 0;    // current event index
  size_t evgen = 0;    // number of generated events
  size_t ngamma = 0;   // number of simulated photons on target
  double tot_xsec = 0; // estimated total cross section
  double flux = 0;     // generated fraction of the bremsstrahlung spectrum
  gen::jpsi_event gen;
};

// initialize the tree for J/Psi event data
void init_tree(TTree& t, mc_event& evbuf) {
  t.Branch("event", &evbuf.event, "event/l");
  t.Branch("evgen", &evbuf.evgen, "evgen/l");
  t.Branch("ngamma", &evbuf.ngamma, "ngamma/l");
  t.Branch("tot_xsec", &evbuf.tot_xsec);
  t.Branch("flux", &evbuf.flux);
  t.Branch("xsec", &evbuf.gen.xsec);
  t.Branch("weight", &evbuf.gen.weight);
  t.Branch("s", &evbuf.gen.s);
  t.Branch("t", &evbuf.gen.t);
  t.Branch("tmin", &evbuf.gen.tmin);
  t.Branch("tmax", &evbuf.gen.tmax);
  t.Branch("W", &evbuf.gen.W);
  t.Branch("beam", &evbuf.gen.beam);
  t.Branch("target", &evbuf.gen.target);
  t.Branch("recoil", &evbuf.gen.recoil);
  t.Branch("jpsi", &evbuf.gen.jpsi);
  t.Branch("positron", &evbuf.gen.positron);
  t.Branch("electron", &evbuf.gen.electron);
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
  mc_event evbuf;
  init_tree(*t, evbuf);

  // init the RNG, use the run number as random seed
  std::shared_ptr<TRandom> rng {std::make_shared<TRandom3>()};
  rng->SetSeed(conf.get<int>("run"));

  // init the photon beam, add some diagnostic histos
  gen::bremsstrahlung photon_gen(settings, "mc/photon_gen", rng);
#if 0
  photon_gen.add_histo(ofile, "Egamma", "Photon Energy",
                       {"Egamma",
                        [](const gen::photon_beam& b) { return b.energy; }, 100,
                        photon_gen.range()});
#endif

  // init the J/Psi generator, add some diagnostic histos
  gen::jpsi jpsi_gen(settings, "mc/jpsi_gen", rng);
  // progress meter
  progress_meter progress{events};

  // generate our events
  double cum_xsec = 0; // cumulative xsec
  while (evbuf.event < events) {
    // generate a photon
    auto photon = photon_gen.generate();
    evbuf.ngamma += 1;
    // generate a new physics event using this photon
    evbuf.gen = jpsi_gen.generate(photon);
    // make sure the physics event is good
    if (!evbuf.gen.good) {
      continue;
    }
    // we have a good generated event, update the generator level counters
    evbuf.evgen += 1;
    cum_xsec += evbuf.gen.xsec;
    evbuf.tot_xsec = cum_xsec / evbuf.evgen;
    // check acceptance
    // TODO
    // store this event
    evbuf.event += 1;
    evbuf.flux = photon.flux;
    t->Fill();
    progress.update();
  }

  t->AutoSave();

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
