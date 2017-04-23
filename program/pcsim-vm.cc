#include <TFile.h>
#include <TRandom3.h>
#include <memory>
#include <pcsim/core/assert.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/event.hh>
#include <pcsim/core/event_out.hh>
#include <pcsim/core/exception.hh>
#include <pcsim/core/framework.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/core/progress_meter.hh>
#include <pcsim/gen/beam/beam.hh>
#include <pcsim/gen/beam/photon.hh>
#include <pcsim/gen/decay/gamma_p_event.hh>
#include <pcsim/gen/process/gamma_p_1X.hh>
#include <pcsim/gen/process/gamma_p_2vmX.hh>
#include <pcsim/physics/decay.hh>

using namespace pcsim;

// util function
std::string to_string_exp(double d) {
  std::stringstream ss;
  ss << std::scientific << d;
  return ss.str();
}

// TODO add proper support/constructor for Pc decay
// TODO move hard scattering info to a gamma_p_event super-class
class vm_event : public event {
public:

  // constructor
  vm_event() = default;
  vm_event(const vm_event&) = default;
  vm_event& operator=(const vm_event&) = default;
  vm_event(const size_t evgen, const double xs, const double w = 1.)
      : event{evgen, xs, w} {}

  int add_photon(const beam::photon_data& photon) {
    photon_index_ = add_daughter(photon.beam(), beam_index());
    scat_index_ = add_daughter(photon.beam(), beam_index());

    W_ = sqrt(photon.W2());
    Q2_ = photon.Q2();
    nu_ = photon.nu();
    x_ = photon.x();
    y_ = photon.y();
    update_cross_section(photon.cross_section());
    update_epsilon(photon.epsilon());

    return photon_index_;
  }
  int add_vm(const process::gamma_p_2vmX_data& vm) {
    vm_index_ = add_daughter(vm.vm(), photon_index(), target_index());
    recoil_index_ = add_daughter(vm.X(), photon_index(), target_index());
    
    t_ = vm.t();
    xv_ = vm.xv();
    Q2plusMv2_ = vm.Q2plusMv2();

    update_cross_section(vm.cross_section());
    update_R(vm.R());

    return vm_index_;
  }

  double W() const {return W_;}
  double W2() const { return W_ * W_; }
  double Q2() const { return Q2_; }
  double nu() const { return nu_; }
  double x() const { return x_; }
  double y() const { return y_; }
  double t() const { return t_; }
  double xv() const { return xv_; }
  double Q2plusMv2() const { return Q2plusMv2_; }
  
  // particle info
  int photon_index() const { return photon_index_; }
  int scat_index() const { return scat_index_; }
  int vm_index() const { return vm_index_; }
  int recoil_index() const { return recoil_index_; }

private:
  double W_ = 0;
  double Q2_ = 0;
  double nu_ = 0;
  double x_ = 0;
  double y_ = 0;

  double t_ = 0;
  double xv_ = 0;
  double Q2plusMv2_ = 0;

  int photon_index_;
  int scat_index_;
  int vm_index_;
  int recoil_index_;
};
  
class vm_event_out : public event_out {
public:
  vm_event_out(std::shared_ptr<TFile> f, const std::string& name)
      : event_out{f, name} {
    create_branches();
  }

  void push(const vm_event& e) {
    W_ = static_cast<float>(e.W());
    Q2_ = static_cast<float>(e.Q2());
    nu_ = static_cast<float>(e.nu());
    x_ = static_cast<float>(e.x());
    y_ = static_cast<float>(e.y());
    t_ = static_cast<float>(e.t());
    Q2plusMv2_ = static_cast<float>(e.Q2plusMv2());
    photon_index_ = static_cast<int16_t>(e.photon_index());
    scat_index_ = static_cast<int16_t>(e.scat_index());
    vm_index_ = static_cast<int16_t>(e.vm_index());
    recoil_index_ = static_cast<int16_t>(e.recoil_index());

    event_out::push(e);
  }

private:
  void create_branches() {
    tree()->Branch("W", &W_);
    tree()->Branch("Q2", &Q2_);
    tree()->Branch("nu", &nu_);
    tree()->Branch("x", &x_);
    tree()->Branch("y", &y_);
    tree()->Branch("t", &t_);
    tree()->Branch("xv", &xv_);
    tree()->Branch("Q2plusMv2", &Q2plusMv2_);
    tree()->Branch("photon_index", &photon_index_);
    tree()->Branch("scat_index", &scat_index_);
    tree()->Branch("vm_index", &vm_index_);
    tree()->Branch("recoil_index", &recoil_index_);
  }

  float W_;
  float Q2_;
  float nu_;
  float x_;
  float y_;
  float t_;
  float xv_;
  float Q2plusMv2_;
  int16_t photon_index_;
  int16_t scat_index_;
  int16_t vm_index_;
  int16_t recoil_index_;
};

class vm_generator : public event_generator<vm_event> {
public:
  using parent_type = event_generator<vm_event>;

  vm_generator(const configuration& cf, const string_path& path,
               std::shared_ptr<TRandom> r)
      : parent_type{cf, path, r}
      , electron_gen{FACTORY_CREATE(beam::primary, conf(), "beam", r)}
      , proton_gen{FACTORY_CREATE(beam::primary, conf(), "target", r)}
      , photon_gen{FACTORY_CREATE(beam::photon, conf(), "photon", r)}
      , vm_gen{FACTORY_CREATE(process::gamma_p_2vmX, conf(), "gamma_p_2vmX", r)}
      , decay_gen{std::make_unique<decay::gamma_p_event>(r)} {
    add(*photon_gen);
    add(*vm_gen);
    }

  protected:
    virtual void generate_event(event_type & e) const {
      beam::data electron = electron_gen->generate();
      beam::data proton = proton_gen->generate();
      beam::photon_data photon = photon_gen->generate(electron, proton);
      if (photon.cross_section() == 0) {
        e.update_cross_section(0);
        return;
      }
      process::gamma_p_2vmX_data vm = vm_gen->generate(photon, proton);
      if (vm.cross_section() == 0) {
        e.update_cross_section(0);
        return;
      }
      e.add_beam(electron.beam());
      e.add_target(proton.beam());
      e.add_photon(photon);
      e.add_vm(vm);
    }

    virtual void build_event(event_type & e) const {
      decay_gen->decay(e);
      e.update_evgen(n_events());
    }

  private:
    std::shared_ptr<beam::primary> electron_gen;
    std::shared_ptr<beam::primary> proton_gen;
    std::shared_ptr<beam::photon> photon_gen;
    std::shared_ptr<process::gamma_p_2vmX> vm_gen;
    std::shared_ptr<decay::gamma_p_event> decay_gen;
};

int run_mc(const configuration& cf, const std::string& output) {
  LOG_INFO("pcsim-vm", "initializing PCSIM-vm");

  // get RNG
  std::shared_ptr<TRandom> r {std::make_shared<TRandom3>()};
  r->SetSeed(cf.get<int>("run"));

  // make output file and buffer
  std::shared_ptr<TFile> ofile{
      std::make_shared<TFile>((output + ".root").c_str(), "recreate")};
  vm_event_out evbuf{ofile, "vm_event"};

  // get event generator
  vm_generator gen{cf, "generator", r};

  // number of requested events:
  const size_t events = cf.get<size_t>("events");

  progress_meter progress{events};

  // loop over events
  while (gen.n_events() < events) {
    evbuf.push(gen.generate());
    progress.update();
  }
  LOG_INFO("pcsim", "Event generation complete");
  LOG_INFO("pcsim", "Total number of generated events: " +
                        std::to_string(gen.n_events()));
  LOG_INFO("pcsim",
           "Total cross section [nb]: " + to_string_exp(gen.cross_section()));

  return 0;
}

MAKE_PCSIM_FRAMEWORK(run_mc)
