#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>
#include <iomanip>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/exception.hh>
#include <pcsim/core/framework.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/core/progress_meter.hh>
#include <pcsim/gen/bremsstrahlung.hh>
#include <pcsim/gen/jpsi.hh>
#include <pcsim/gen/pc.hh>
#include <pcsim/gen/spectrometer.hh>
#include <sstream>

using namespace pcsim;

// util function
std::string to_string_exp(double d) {
  std::stringstream ss;
  ss << std::scientific << d;
  return ss.str();
}

// notes:
//  * For a give number of incoming electrons Ne, the number of good
//    bremsstrahlung photons Nb that for the process of interest is given by: 
//      Nb = Ne * flux * (evgen/ngamma)
//  * When binning using the cross section as event weight w_i, the differential
//    cross section in a bin is given by: 
//      dsigma = sum_i(w_i) / evgen / delta
//    with delta the bin width
struct mc_event {
  uint64_t event = 0;  // current event index
  uint64_t evgen = 0;  // number of generated events
  uint64_t ngamma = 0; // number of simulated photons on target
  uint64_t nHMS = 0;   // number of accepted HMS tracks
  uint64_t nSHMS = 0;  // number of accepted SHMS tracks
  double xsec_gen = 0; // estimated total cross section for all generated events
  double xsec_acc = 0; // estimated total cross section for accepted evens
  double flux = 0;     // generated fraction of the bremsstrahlung spectrum
  gen::jpsi_event gen;
  gen::spec_track HMS;
  gen::spec_track SHMS;
};

struct mc_controller {
  const configuration conf;
  const int run;
  const size_t events;
  const bool spectrometer;
  const bool pc_mode;
  std::shared_ptr<TFile> ofile;
  TTree* tree;
  std::shared_ptr<TRandom> rng;
  progress_meter progress;
  mc_event ev;
  double cum_xsec_gen; // cumulative generated cross section
  double cum_xsec_acc; // cumulative accepted cross section

  mc_controller(const ptree& settings, const std::string& output)
      : conf{settings, "mc"}
      , run{conf.get<int>("run")}
      , events{conf.get<size_t>("events")}
      , spectrometer{conf.get<bool>("spectrometer")}
      , pc_mode{conf.get<bool>("pc_mode")}
      , ofile{std::make_shared<TFile>((output + ".root").c_str(), "recreate")}
      , tree{new TTree{"jpsi_event", "J/Psi Event Data"}}
      , rng{std::make_shared<TRandom3>()}
      , progress{events}
      , cum_xsec_gen{0}
      , cum_xsec_acc{0} {
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
    tree->Branch("nHMS", &ev.nHMS, "nHMS/l");
    tree->Branch("nSHMS", &ev.nSHMS, "nSHMS/l");
    tree->Branch("xsec_gen", &ev.xsec_gen);
    tree->Branch("xsec_acc", &ev.xsec_gen);
    tree->Branch("flux", &ev.flux);
    tree->Branch("xsec", &ev.gen.xsec);
    tree->Branch("weight", &ev.gen.weight);
    tree->Branch("branching", &ev.gen.branching);
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
    // spectrometer data
    tree->Branch("HMS.track", &ev.HMS.track);
    tree->Branch("HMS.charge", &ev.HMS.charge);
    tree->Branch("HMS.p", &ev.HMS.p);
    tree->Branch("HMS.thx", &ev.HMS.thx);
    tree->Branch("HMS.thy", &ev.HMS.thy);
    tree->Branch("SHMS.track", &ev.SHMS.track);
    tree->Branch("SHMS.charge", &ev.SHMS.charge);
    tree->Branch("SHMS.p", &ev.SHMS.p);
    tree->Branch("SHMS.thx", &ev.SHMS.thx);
    tree->Branch("SHMS.thy", &ev.SHMS.thy);
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
    cum_xsec_gen += ev.gen.xsec;
    ev.xsec_gen = cum_xsec_gen / ev.evgen;
  }

  void book_event() {
    ev.event += 1;
    // update the estimated accepted cross section
    cum_xsec_acc += ev.gen.xsec;
    ev.xsec_acc = cum_xsec_acc / ev.evgen;
    // book the data
    tree->Fill();
    progress.update();
  }
};

int run_mc(const ptree& settings, const std::string& output) { 
  LOG_INFO("MC", "initializing PCSIM"); 

  // mc settings
  mc_controller mc{settings, output};

  // init the photon beam, add some diagnostic histos
  gen::bremsstrahlung photon_gen(settings, "mc/photon_gen", mc.rng);

  // init the J/Psi generator
  gen::jpsi jpsi_gen{settings, "mc/jpsi_gen", mc.rng};

  // init the Pc generator
  gen::pc pc_gen{settings, "mc/pc_gen", mc.rng};

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
    if (!mc.pc_mode) {
      mc.ev.gen = jpsi_gen.generate(photon);
    } else {
      mc.ev.gen = pc_gen.generate(photon);
    }

    // make sure the physics event is good
    if (!mc.ev.gen.good) {
      continue;
    }
    // we have a good generated event, update the generator level counters
    mc.record_xsec();

    // check acceptance
    if (mc.spectrometer) {
      mc.ev.HMS = HMS.check(mc.ev.gen.positron, 1);
      if (mc.ev.HMS.accept) {
        mc.ev.nHMS++;
      };
      mc.ev.SHMS = SHMS.check(mc.ev.gen.electron, -1);
      if (mc.ev.SHMS.accept) {
        mc.ev.nSHMS++;
      };
      if (!mc.ev.HMS.accept || !mc.ev.SHMS.accept) {
        // not in acceptance, try again
        continue;
      }
    }
    // store this event
    mc.book_event();
  }

  // give some stats
  LOG_INFO("MC", "Event generation complete");
  LOG_INFO("MC", "Total number of photons on target: " + std::to_string(mc.ev.ngamma));
  LOG_INFO("MC",
           "Total number of generated events: " + std::to_string(mc.ev.evgen));
  LOG_INFO("MC",
           "Total cross section [nb]: " + to_string_exp(mc.ev.xsec_gen));
  if (mc.spectrometer) {
    LOG_INFO("MC", "Total number accepted positrons in the HMS: " +
                       std::to_string(mc.ev.nHMS));
    LOG_INFO("MC", "Total number accepted electrons in the SHMS: " +
                       std::to_string(mc.ev.nSHMS));
    LOG_INFO("MC", "Number of accepted events: " + std::to_string(mc.ev.event));
    LOG_INFO("MC", "Accepted cross section [nb]: " + to_string_exp(mc.ev.xsec_acc));
  }

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
