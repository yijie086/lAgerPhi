#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>
#include <fstream>
#include <iomanip>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/exception.hh>
#include <pcsim/core/framework.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/core/progress_meter.hh>
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

struct mc_event {
  uint64_t event = 0;   // current event index
  uint64_t evgen = 0;   // number of generated events
  uint64_t ntrials = 0; // number of sampled phase space points
  uint64_t nHMS = 0;    // number of accepted HMS tracks
  uint64_t nSHMS = 0;   // number of accepted SHMS tracks
  double xsec_gen = 0;  // total cross section for all generated events
  double xsec_acc = 0;  // total cross section for accepted evens
  gen::jpsi_event gen;
  gen::spec_track HMS;
  gen::spec_track SHMS;
};

struct mc_controller {
  const configuration conf;
  const int run;
  const size_t events;
  const std::string generator;
  const std::string acceptance;
  std::ofstream logfile;
  std::shared_ptr<TFile> ofile;
  TTree* tree;
  std::shared_ptr<TRandom> rng;
  progress_meter progress;
  mc_event ev;
  double volume; // sum of the phase space volume for all events

  mc_controller(const configuration& cf, const std::string& output)
      : conf{cf}
      , run{conf["run"]}
      , events{conf["events"]}
      , generator{conf.get<std::string>("generator/type")}
      , acceptance{conf.get<std::string>("acceptance/type")}
      , logfile(output + ".log")
      , ofile{std::make_shared<TFile>((output + ".root").c_str(), "recreate")}
      , tree{new TTree{"jpsi_event", "J/Psi Event Data"}}
      , rng{std::make_shared<TRandom3>()}
      , progress{events}
      , volume{0} {
    // redirect logger to use the log file
    global::logger.set_output(logfile);
    ofile->cd();
    init_tree();
    // use the run number as random seed for the RNG
    rng->SetSeed(run);
  }
  ~mc_controller() {
    tree->AutoSave();
    // reset the logger to use standard output
    global::logger.set_output(std::cout);
  }
  // initialize the tree for J/Psi event data
  void init_tree() {
    tree->Branch("event", &ev.event, "event/l");
    tree->Branch("evgen", &ev.evgen, "evgen/l");
    tree->Branch("ntrials", &ev.evgen, "ntrials/l");
    tree->Branch("nHMS", &ev.nHMS, "nHMS/l");
    tree->Branch("nSHMS", &ev.nSHMS, "nSHMS/l");
    tree->Branch("xsec_gen", &ev.xsec_gen);
    tree->Branch("xsec_acc", &ev.xsec_acc);
    tree->Branch("weight", &ev.gen.weight);
    tree->Branch("s", &ev.gen.s);
    tree->Branch("W", &ev.gen.W);
    tree->Branch("t", &ev.gen.t);
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
  void record_geninfo() {
    // update the counters
    ev.evgen += 1;
    ev.ntrials += ev.gen.ntrials;
    volume += ev.gen.volume;
    // update the estimated total cross section
    ev.xsec_gen = volume / ev.ntrials;
  }

  void book_event() {
    ev.event += 1;
    // update the estimated accepted cross section
    ev.xsec_acc = volume * ev.event / ev.evgen / ev.ntrials;
    // book the data
    tree->Fill();
    progress.update();
  }
};

int run_mc(const configuration& cf, const std::string& output) {
  LOG_INFO("pcsim", "initializing PCSIM");

  // mc settings
  mc_controller mc{cf, output};

  // init the physics generator

  std::unique_ptr<gen::pc> pc_gen;
  std::unique_ptr<gen::jpsi> jpsi_gen;

  if (mc.generator == "pc") {
    pc_gen = std::make_unique<gen::pc>(mc.conf, "generator", mc.rng);
  } else if (mc.generator == "jpsi") {
    jpsi_gen = std::make_unique<gen::jpsi>(mc.conf, "generator", mc.rng);
  } else {
    throw mc.conf.value_error("generator/type");
  }

  // init the spectrometers
  std::unique_ptr<gen::spectrometer> HMS;
  std::unique_ptr<gen::spectrometer> SHMS;
  if (mc.acceptance == "2arm") {
    HMS =
        std::make_unique<gen::spectrometer>(mc.conf, "acceptance/HMS", mc.rng);
    SHMS =
        std::make_unique<gen::spectrometer>(mc.conf, "acceptance/SHMS", mc.rng);
  }

  // generate our events
  while (mc.ev.event < mc.events) {

    // generate a new physics event using this photon
    if (pc_gen) {
      mc.ev.gen = pc_gen->generate();
    } else {
      mc.ev.gen = jpsi_gen->generate();
    }
    // we have a good generated event, update the generator level counters
    mc.record_geninfo();

    // check acceptance
    if (HMS && SHMS) {
      mc.ev.HMS = HMS->check(mc.ev.gen.electron, -1);
      if (mc.ev.HMS.accept) {
        mc.ev.nHMS++;
      };
      mc.ev.SHMS = SHMS->check(mc.ev.gen.positron, 1);
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
  LOG_INFO("pcsim", "Event generation complete");
  LOG_INFO("pcsim",
           "Total number of generated events: " + std::to_string(mc.ev.evgen));
  LOG_INFO("pcsim",
           "Total cross section [nb]: " + to_string_exp(mc.ev.xsec_gen));
  if (mc.acceptance != "4pi") {
    LOG_INFO("pcsim", "Total number accepted positrons in the HMS: " +
                          std::to_string(mc.ev.nHMS));
    LOG_INFO("pcsim", "Total number accepted electrons in the SHMS: " +
                          std::to_string(mc.ev.nSHMS));
    LOG_INFO("pcsim",
             "Number of accepted events: " + std::to_string(mc.ev.event));
    LOG_INFO("pcsim",
             "Accepted cross section [nb]: " + to_string_exp(mc.ev.xsec_acc));
    LOG_INFO("pcsim", "Acceptance [%]: " +
                          std::to_string(static_cast<float>(mc.ev.event) /
                                         mc.ev.evgen * 100.f));
  }

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
