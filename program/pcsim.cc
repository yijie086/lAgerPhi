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
#include <pcsim/gen/spectrometer.hh>

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
  uint64_t event = 0;    // current event index
  uint64_t evgen = 0;    // number of generated events
  uint64_t ngamma = 0;   // number of simulated photons on target
  double tot_xsec = 0; // estimated total cross section
  double flux = 0;     // generated fraction of the bremsstrahlung spectrum
  gen::jpsi_event gen;
};

struct mc_controller {
  const configuration conf;
  const int run;
  const size_t events;
  const bool spectrometer;
  std::shared_ptr<TFile> ofile;
  TTree* tree;
  std::shared_ptr<TRandom> rng;
  progress_meter progress;
  mc_event ev;
  double cum_xsec; // cumulative cross section

  mc_controller(const ptree& settings, const std::string& output)
      : conf{settings, "mc"}
      , run{conf.get<int>("run")}
      , events{conf.get<size_t>("events")}
      , spectrometer{conf.get<bool>("spectrometer")}
      , ofile{std::make_shared<TFile>((output + ".root").c_str(), "recreate")}
      , tree{new TTree{"jpsi_event", "J/Psi Event Data"}}
      , rng{std::make_shared<TRandom3>()}
      , progress{events}
      , cum_xsec{0} {
    ofile->cd();
    init_tree();
    // use the run number as random seed for the RNG
    rng->SetSeed(run);
  }
  ~mc_controller() { tree->AutoSave(); }
  // initialize the tree for J/Psi event data
  void init_tree() {
    tree->Branch("event", &ev.event, "event/l");
    tree->Branch("evgen", &ev.evgen, "evgen/l");
    tree->Branch("ngamma", &ev.ngamma, "ngamma/l");
    tree->Branch("tot_xsec", &ev.tot_xsec);
    tree->Branch("flux", &ev.flux);
    tree->Branch("xsec", &ev.gen.xsec);
    tree->Branch("weight", &ev.gen.weight);
    tree->Branch("s", &ev.gen.s);
    tree->Branch("t", &ev.gen.t);
    tree->Branch("tmin", &ev.gen.tmin);
    tree->Branch("tmax", &ev.gen.tmax);
    tree->Branch("W", &ev.gen.W);
    tree->Branch("beam", &ev.gen.beam);
    tree->Branch("target", &ev.gen.target);
    tree->Branch("recoil", &ev.gen.recoil);
    tree->Branch("jpsi", &ev.gen.jpsi);
    tree->Branch("positron", &ev.gen.positron);
    tree->Branch("electron", &ev.gen.electron);
  }
  void record_photon(const gen::photon_beam& photon) {
    // update the counter
    ev.ngamma += 1;
    // store the flux
    ev.flux = photon.flux;
  }
  void record_xsec() {
    // update the counter
    ev.evgen += 1;
    // update the estimated total cross section
    cum_xsec += ev.gen.xsec;
    ev.tot_xsec = cum_xsec / ev.evgen;
  }

  void book_event() {
    ev.event += 1;
    tree->Fill();
    progress.update();
  }
};

int run_mc(const ptree& settings, const std::string& output) { 
  LOG_INFO("main", "PCSIM"); 

  // mc settings
  mc_controller mc{settings, output};

  // init the photon beam, add some diagnostic histos
  gen::bremsstrahlung photon_gen(settings, "mc/photon_gen", mc.rng);

  // init the J/Psi generator
  gen::jpsi jpsi_gen{settings, "mc/jpsi_gen", mc.rng};

  // init the spectrometers
  gen::spectrometer HMS{settings, "mc/HMS", mc.rng};
  gen::spectrometer SHMS{settings, "mc/SHMS", mc.rng};

  // generate our events
  double cum_xsec = 0; // cumulative xsec
  while (mc.ev.event < mc.events) {

    // generate a photon
    auto photon = photon_gen.generate();
    mc.record_photon(photon);

    // generate a new physics event using this photon
    mc.ev.gen = jpsi_gen.generate(photon);

    // make sure the physics event is good
    if (!mc.ev.gen.good) {
      continue;
    }
    // we have a good generated event, update the generator level counters
    mc.record_xsec();

    // check acceptance
    if (mc.spectrometer) {
      if (!HMS.check(mc.ev.gen.positron).accept ||
          !SHMS.check(mc.ev.gen.electron).accept) {
        continue;
      }
    }
    // store this event
    mc.book_event();
  }

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
