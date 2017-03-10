#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/event.hh>
#include <pcsim/core/exception.hh>
#include <pcsim/core/framework.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/core/progress_meter.hh>
#include <pcsim/gen/beam.hh>
#include <pcsim/gen/exclusive.hh>
#include <pcsim/gen/photon.hh>

using namespace pcsim;

// util function
std::string to_string_exp(double d) {
  std::stringstream ss;
  ss << std::scientific << d;
  return ss.str();
}
struct vm_event : event {
  int evgen;
  double W = 0;
  double Q2 = 0;
  double nu = 0;
  double x = 0;
  double y = 0;
  double t = 0;
  double xv = 0;
  double Q2plusMv2 = 0;
  TLorentzVector beam;
  TLorentzVector target;
  TLorentzVector scat;
  TLorentzVector photon;
  TLorentzVector recoil;
  TLorentzVector vm;
  TLorentzVector positron;
  TLorentzVector electron;

  int scat_index = -1;
  int photon_index = -1;
  int vm_index = -1;
  int recoil_index = -1;
  int positron_index = -1;
  int electron_index = -1;

  vm_event() { part.reserve(30); }

  void link(TFile& f) {
    f.cd();
    t_ = new TTree("vm_event", "Upsilon Event Data");
    t_->Branch("event", &event_);
    t_->Branch("evgen", &evgen);
    t_->Branch("cross_section", &cross_section);
    t_->Branch("weight", &weight);
    t_->Branch("W", &W);
    t_->Branch("Q2", &Q2);
    t_->Branch("nu", &nu);
    t_->Branch("x", &x);
    t_->Branch("y", &y);
    t_->Branch("t", &t);
    t_->Branch("xv", &xv);
    t_->Branch("Q2plusMv2", &Q2plusMv2);
    t_->Branch("beam", &beam);
    t_->Branch("target", &target);
    t_->Branch("scat", &scat);
    t_->Branch("photon", &photon);
    t_->Branch("vm", &vm);
    t_->Branch("recoil", &recoil);
    t_->Branch("positron", &positron);
    t_->Branch("electron", &electron);
  }
  void write() {
    if (t_) {
      t_->Fill();
      // increment event index;
      event_ += 1;
    }
  }

  vm_event& operator=(const vm_event& rhs) {
    static_cast<event&>(*this) = rhs;
    evgen = rhs.evgen;
    W = rhs.W;
    Q2 = rhs.Q2;
    nu = rhs.nu;
    x = rhs.x;
    y = rhs.y;
    t = rhs.t;
    xv = rhs.xv;
    Q2plusMv2 = rhs.Q2plusMv2;
    scat_index = rhs.scat_index;
    photon_index = rhs.photon_index;
    vm_index = rhs.vm_index;
    recoil_index = rhs.recoil_index;
    positron_index = rhs.positron_index;
    electron_index = rhs.electron_index;

    beam = part[beam_index].mom;
    target = part[target_index].mom;
    scat = part[scat_index].mom;
    photon = part[photon_index].mom;
    vm = part[vm_index].mom;
    recoil = part[recoil_index].mom;
    electron = part[electron_index].mom;
    positron = part[positron_index].mom;
    return *this;
  }

  ~vm_event() {
    if (t_) {
      t_->AutoSave();
    }
  }

private:
  int event_ = 0;
  TTree* t_ = 0;
};

class vm_generator : public event_generator<vm_event> {
public:
  using parent_type = event_generator<vm_event>;

  vm_generator(const configuration& cf, const string_path& path,
                    std::shared_ptr<TRandom> r)
      : parent_type{cf, path, r}
      , electron_gen{std::make_unique<gen::beam>(conf(), "beam", r)}
      , proton_gen{std::make_unique<gen::beam>(conf(), "target", r)}
      , photon_gen{conf().get<std::string>("photon/type") == "vphoton"
                       ? (gen::photon*)new gen::vphoton{conf(), "photon", r}
                       : (gen::photon*)new gen::bremsstrahlung{conf(), "photon",
                                                               r}}
      , excl_gen{std::make_unique<gen::brodsky_tchannel>(conf(), "excl", r)} {
    add_max_cross_section(photon_gen->max_cross_section());
    add_max_cross_section(excl_gen->max_cross_section());
  }

protected:
  virtual void generate_event(event_type& e) const {
    particle e_in = electron_gen->generate();
    particle p_in = proton_gen->generate();

    // create lepton/proton beam
    e.part.push_back({e_in, mc_particle::INITIAL});
    e.part.push_back({p_in, mc_particle::INITIAL});
    e.beam_index = 0;
    e.target_index = 1;

    // create virtual photon and scattered lepton
    gen::photon_data photon_data = photon_gen->generate(e_in, p_in);
    e.W = sqrt(photon_data.W2);
    e.Q2 = photon_data.Q2;
    e.nu = photon_data.nu;
    e.x = photon_data.x;
    e.y = photon_data.y;
    e.part.push_back({photon_data.scat, mc_particle::FINAL});
    e.part.push_back({photon_data.photon, mc_particle::VIRTUAL});
    e.scat_index = 2;
    e.photon_index = 3;
    e.cross_section = photon_data.cross_section;
    e.phase_space = photon_data.phase_space;
    if (e.cross_section <= 0) {
      return;
    }

    // create vm and recoil
    gen::exclusive_data excl_data = excl_gen->generate(photon_data, p_in);
    e.t = excl_data.t;
    e.xv = excl_data.xv;
    e.Q2plusMv2 = excl_data.Q2plusMv2;
    e.part.push_back({excl_data.vm, mc_particle::VIRTUAL});
    e.part.push_back({excl_data.recoil, mc_particle::FINAL});
    e.vm_index = 4;
    e.recoil_index = 5;
    e.cross_section *= excl_data.cross_section;
    e.phase_space *= excl_data.phase_space;
    if (e.cross_section <= 0) {
      return;
    }

    e.evgen = nevents();
  }

  virtual void build_event(event_type& e) const {
    ; // no-op
  }

private:
  std::unique_ptr<gen::beam> electron_gen;
  std::unique_ptr<gen::beam> proton_gen;
  std::unique_ptr<gen::photon> photon_gen;
  std::unique_ptr<gen::exclusive> excl_gen;
};

int run_mc(const configuration& cf, const std::string& output) {
  LOG_INFO("pcsim-vm", "initializing PCSIM-vm");

  // get RNG
  std::shared_ptr<TRandom> r {std::make_shared<TRandom3>()};
  r->SetSeed(cf.get<int>("run"));

  // make output file and buffer
  TFile ofile{(output + ".root").c_str(), "recreate"};
  vm_event evbuf;
  evbuf.link(ofile);

  // get event generator
  vm_generator gen{cf, "generator", r};

  // number of requested events:
  const size_t events = cf.get<size_t>("events");

  progress_meter progress{events};

  // loop over events
  while (gen.nevents() < events) {
    evbuf = gen.generate();
    evbuf.write();
    progress.update();
  }
  LOG_INFO("pcsim", "Event generation complete");
  LOG_INFO("pcsim",
           "Total number of generated events: " +
               std::to_string(gen.nevents()));
  LOG_INFO("pcsim",
           "Total cross section [nb]: " + to_string_exp(gen.cross_section()));

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
