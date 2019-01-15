#ifndef PCSIM_GEN_LP_GAMMA_EVENT_LOADED
#define PCSIM_GEN_LP_GAMMA_EVENT_LOADED

#include <fstream>
#include <pcsim/core/event.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/gen/beam/photon.hh>
#include <pcsim/gen/beam/primary.hh>

namespace pcsim {

// =============================================================================
// Initial reaction data for a lepton-proton reaction mediated by a virtual
// photon.
// =============================================================================
class lp_gamma_data : public generator_data {
public:
  lp_gamma_data(const double xs = 0) : generator_data{xs} {}
  lp_gamma_data(const beam::primary& l, const beam::primary& p,
                const beam::photon& gamma);
  lp_gamma_data(const lp_gamma_data&) = default;
  lp_gamma_data& operator=(const lp_gamma_data&) = default;

  beam::primary& lepton() { return lepton_; }
  const beam::primary& lepton() const { return lepton_; }
  beam::primary& target() { return target_; }
  const beam::primary& target() const { return target_; }
  beam::photon& photon() { return photon_; }
  const beam::photon& photon() const { return photon_; }

private:
  beam::primary lepton_;
  beam::primary target_;
  beam::photon photon_;
};

// =============================================================================
// full lp_gamma event
// =============================================================================
class lp_gamma_event : public event {
public:
  lp_gamma_event(const lp_gamma_event&) = default;
  lp_gamma_event& operator=(const lp_gamma_event&) = default;
  explicit lp_gamma_event(const double xs = 1., const double w = 1.,
                          const double R = 0., const double epsilon = 0.);
  lp_gamma_event(const lp_gamma_data& initial, const double xs,
                 const double w = 1., const double R = 0.);

  // ===========================================================================
  // PHOTON INFO
  double W() const { return W_; }
  double W2() const { return W_ * W_; }
  double Q2() const { return Q2_; }
  double nu() const { return nu_; }
  double x() const { return x_; }
  double y() const { return y_; }
  double epsilon() const { return epsilon_; }
  double R() const { return R_; }
  double epsilon_R() const { return epsilon_ * R_; }

  // modifiers
  void update_R(const double R) { R_ = R; }
  void update_epsilon(const double epsilon) { epsilon_ = epsilon; }

  // ===========================================================================
  // PARTICLE INFO
  //
  // add photon and scat
  int add_photon(const beam::photon& gamma);
  // add leading particle
  int add_leading(const particle& part, const int parent1,
                  const int parent2 = -1);
  int add_leading(const particle& part);
  // add recoil particle
  int add_recoil(const particle& part, const int parent1,
                 const int parent2 = -1);
  int add_recoil(const particle& part);

  // get photon, scat, leading and recoil info
  particle& photon();
  const particle& photon() const;
  particle& scat();
  const particle& scat() const;
  particle& leading();
  const particle& leading() const;
  particle& recoil();
  const particle& recoil() const;

  int photon_index() const { return photon_index_; }
  int scat_index() const { return scat_index_; }
  int leading_index() const { return leading_index_; }
  int recoil_index() const { return recoil_index_; }

  // ===========================================================================
  // DETECTED PARTICLE INFO
  //
  // set the relevant detected particle index
  void update_detected_index(const int index);

  // get photon, scat, leading and recoil info
  detected_particle& detected_photon();
  const detected_particle& detected_photon() const;
  detected_particle& detected_scat();
  const detected_particle& detected_scat() const;
  detected_particle& detected_leading();
  const detected_particle& detected_leading() const;
  detected_particle& detected_recoil();
  const detected_particle& detected_recoil() const;

  int detected_photon_index() const { return rc_photon_index_; }
  int detected_scat_index() const { return rc_scat_index_; }
  int detected_leading_index() const { return rc_leading_index_; }
  int detected_recoil_index() const { return rc_recoil_index_; }

  // ===========================================================================
  // leading produced particle kinematics
  //
  // Mandelstam t
  double t() const;
  // modified bjorken x (xv = (Q2 + M2)/(2Mnu))
  double xv() const;
  // modified Q2
  double Q2plusM2() const;

  // ===========================================================================
  // DETECTED KINEMATICS
  //
  double detected_W() const;
  double detected_W2() const;
  double detected_Q2() const;
  double detected_nu() const;
  double detected_x() const;
  double detected_y() const;
  double detected_t() const;
  double detected_xv() const;
  double detected_Q2plusM2() const;

private:
  // photon info
  double W_{0.};
  double Q2_{0.};
  double nu_{0.};
  double x_{0.};
  double y_{0.};
  double epsilon_{0.};
  double R_{0.};
  int photon_index_{-1};
  int scat_index_{-1};

  // leading and recoil produced particle info
  int leading_index_{-1};
  int recoil_index_{-1};

  // and reconstructed
  int rc_photon_index_{-1};
  int rc_scat_index_{-1};
  int rc_leading_index_{-1};
  int rc_recoil_index_{-1};
};

// =============================================================================
// lp_gamma_out custom event output
// =============================================================================
class lp_gamma_out : public event_out {
public:
  lp_gamma_out(std::shared_ptr<TFile> f, std::ofstream& olund,
               const std::string& name);
  void push(const lp_gamma_event& e);
  void push(const std::vector<lp_gamma_event>& e);

private:
  void create_branches();

  // photon kinematics
  float W_;
  float Q2_;
  float nu_;
  float x_;
  float y_;
  float epsilon_;
  float R_;
  // leading particle kinematics
  float t_;
  float xv_;
  float Q2plusM2_;
  // particle indices
  int16_t photon_index_;
  int16_t scat_index_;
  int16_t leading_index_;
  int16_t recoil_index_;
  // and reconstructed
  float rc_W_;
  float rc_Q2_;
  float rc_nu_;
  float rc_x_;
  float rc_y_;
  float rc_t_;
  float rc_xv_;
  float rc_Q2plusM2_;
  int16_t rc_photon_index_;
  int16_t rc_scat_index_;
  int16_t rc_leading_index_;
  int16_t rc_recoil_index_;
};

// =============================================================================
// LP_GAMMA_DATA IMPLEMENTATION
// =============================================================================
inline lp_gamma_data::lp_gamma_data(const beam::primary& l,
                                    const beam::primary& p,
                                    const beam::photon& gamma)
    : generator_data{gamma.cross_section()}
    , lepton_{l}
    , target_{p}
    , photon_{gamma} {}

// =============================================================================
// LP_GAMMA_EVENT IMPLEMENTATION
// =============================================================================

// constructors
inline lp_gamma_event::lp_gamma_event(const double xs, const double w,
                                      const double R, const double epsilon)
    : event{xs, w}, R_{R}, epsilon_{epsilon} {}
inline lp_gamma_event::lp_gamma_event(const lp_gamma_data& initial,
                                      const double xs, const double w,
                                      const double R)
    : event{xs, w}, R_{R}, epsilon_{initial.photon().epsilon()} {
  add_beam(initial.lepton().beam());
  add_target(initial.target().beam());
  // thes step also applies the photon cross section!
  add_photon(initial.photon());
}

// set particles
inline int lp_gamma_event::add_photon(const beam::photon& gamma) {
  photon_index_ = add_daughter(gamma.beam(), beam_index());
  scat_index_ = add_daughter(gamma.scat(), beam_index());

  W_ = sqrt(gamma.W2());
  Q2_ = gamma.Q2();
  nu_ = gamma.nu();
  x_ = gamma.x();
  y_ = gamma.y();

  update_cross_section(gamma.cross_section());

  return photon_index();
}
inline int lp_gamma_event::add_leading(const particle& part, const int parent1,
                                       const int parent2) {
  leading_index_ = add_daughter(part, parent1, parent2);
  return leading_index();
}
inline int lp_gamma_event::add_leading(const particle& part) {
  return add_leading(part, photon_index(), target_index());
}
inline int lp_gamma_event::add_recoil(const particle& part, const int parent1,
                                      const int parent2) {
  recoil_index_ = add_daughter(part, parent1, parent2);
  return recoil_index();
}
inline int lp_gamma_event::add_recoil(const particle& part) {
  return add_recoil(part, photon_index(), target_index());
}

// get photon, scat, leading and recoil info
inline particle& lp_gamma_event::photon() {
  tassert(photon_index_ >= 0, "trying to access photon data, but no photon "
                              "data present in the event");
  return part(photon_index());
}
inline const particle& lp_gamma_event::photon() const {
  tassert(photon_index_ >= 0, "trying to access photon data, but no photon "
                              "data present in the event");
  return part(photon_index());
}
inline particle& lp_gamma_event::scat() {
  tassert(scat_index_ >= 0, "trying to access scat data, but no scat "
                            "data present in the event");
  return part(scat_index());
}
inline const particle& lp_gamma_event::scat() const {
  tassert(scat_index_ >= 0, "trying to access scat data, but no scat "
                            "data present in the event");
  return part(scat_index());
}
inline particle& lp_gamma_event::leading() {
  tassert(leading_index_ >= 0, "trying to access leading data, but no leading "
                               "data present in the event");
  return part(leading_index());
}
inline const particle& lp_gamma_event::leading() const {
  tassert(leading_index_ >= 0, "trying to access leading data, but no leading "
                               "data present in the event");
  return part(leading_index());
}
inline particle& lp_gamma_event::recoil() {
  tassert(recoil_index_ >= 0, "trying to access recoil data, but no recoil "
                              "data present in the event");
  return part(recoil_index());
}
inline const particle& lp_gamma_event::recoil() const {
  tassert(recoil_index_ >= 0, "trying to access recoil data, but no recoil "
                              "data present in the event");
  return part(recoil_index());
}
// detected particle, also set the relevant detected index
inline void lp_gamma_event::update_detected_index(const int index) {
  const auto& dp = detected(index);
  if (dp.generated().index() == scat_index()) {
    LOG_JUNK2("lp_gamma_event", "Found reconstructed scattered lepton");
    rc_scat_index_ = index;
    // also reconstruct the photon
    // status set to -1 for purely reconstucted particle
    if (beam_index() >= 0 && photon_index() >= 0) {
      LOG_JUNK2("lp_gamma_event", "Adding reconstructed photon data");
      rc_photon_index_ =
          add_detected({photon(), beam().p() - detected_scat().p(), -1});
    } else {
      LOG_JUNK2("lp_gamma_event", "Unable to add reconstructed photon data");
    }
  } else if (dp.generated().index() == photon_index()) {
    LOG_JUNK2("lp_gamma_event", "Found reconstructed scattered photon");
    rc_photon_index_ = index;
  } else if (dp.generated().index() == leading_index()) {
    LOG_JUNK2("lp_gamma_event", "Found reconstructed leading particle");
    rc_leading_index_ = index;
  } else if (dp.generated().index() == recoil_index()) {
    LOG_JUNK2("lp_gamma_event", "Found reconstructed recoil particle");
    rc_recoil_index_ = index;
  }
}

inline detected_particle& lp_gamma_event::detected_photon() {
  tassert(detected_photon_index() >= 0,
          "trying to access detected_photon data, but no detected_photon "
          "data present in the event");
  return detected(detected_photon_index());
}
inline const detected_particle& lp_gamma_event::detected_photon() const {
  tassert(detected_photon_index() >= 0,
          "trying to access detected_photon data, but no detected_photon "
          "data present in the event");
  return detected(detected_photon_index());
}
inline detected_particle& lp_gamma_event::detected_scat() {
  tassert(detected_scat_index() >= 0,
          "trying to access detected_scat data, but no detected_scat "
          "data present in the event");
  return detected(detected_scat_index());
}
inline const detected_particle& lp_gamma_event::detected_scat() const {
  tassert(detected_scat_index() >= 0,
          "trying to access detected_scat data, but no detected_scat "
          "data present in the event");
  return detected(detected_scat_index());
}
inline detected_particle& lp_gamma_event::detected_leading() {
  tassert(detected_leading_index() >= 0,
          "trying to access detected_leading data, but no detected_leading "
          "data present in the event");
  return detected(detected_leading_index());
}
inline const detected_particle& lp_gamma_event::detected_leading() const {
  tassert(detected_leading_index() >= 0,
          "trying to access detected_leading data, but no detected_leading "
          "data present in the event");
  return detected(detected_leading_index());
}
inline detected_particle& lp_gamma_event::detected_recoil() {
  tassert(detected_recoil_index() >= 0,
          "trying to access detected_recoil data, but no detected_recoil "
          "data present in the event");
  return detected(detected_recoil_index());
}
inline const detected_particle& lp_gamma_event::detected_recoil() const {
  tassert(detected_recoil_index() >= 0,
          "trying to access detected_recoil data, but no detected_recoil "
          "data present in the event");
  return detected(detected_recoil_index());
}

// leading produced particle kinematics
inline double lp_gamma_event::t() const {
  if (leading_index() < 0 || photon_index() < 0) {
    return 0.;
  }
  return (leading().p() - photon().p()).M2();
}
inline double lp_gamma_event::xv() const {
  if (leading_index() < 0 || photon_index() < 0 || target_index() < 0) {
    return 0.;
  }
  return (Q2() + leading().mass2()) / (2 * target().mass() * nu());
}
inline double lp_gamma_event::Q2plusM2() const {
  if (leading_index() < 0 || photon_index() < 0) {
    return 0.;
  }
  return Q2() + leading().mass2();
}

// detected kinematics
inline double lp_gamma_event::detected_W() const { return sqrt(detected_W2()); }
inline double lp_gamma_event::detected_W2() const {
  if (detected_photon_index() < 0 || target_index() < 0) {
    return 0.;
  }
  return (detected_photon().p() + target().p()).M2();
}
inline double lp_gamma_event::detected_Q2() const {
  if (detected_photon_index() < 0) {
    return 0.;
  }
  return -detected_photon().p().M2();
}
inline double lp_gamma_event::detected_nu() const {
  if (detected_photon_index() < 0 || target_index() < 0) {
    return 0.;
  }
  return (detected_photon().p()).Dot(target().p()) / target().mass();
}
inline double lp_gamma_event::detected_x() const {
  if (detected_photon_index() < 0 || target_index() < 0) {
    return 0.;
  }
  return detected_Q2() / (2 * target().mass() * detected_nu());
}
inline double lp_gamma_event::detected_y() const {
  if (detected_photon_index() < 0 || target_index() < 0 || beam_index() < 0) {
    return 0.;
  }
  return (detected_photon().p()).Dot(target().p()) /
         (beam().p()).Dot(target().p());
}
inline double lp_gamma_event::detected_t() const {
  if (detected_photon_index() < 0 || detected_leading_index() < 0) {
    return 0.;
  }
  return (detected_leading().p() - detected_photon().p()).M2();
}
inline double lp_gamma_event::detected_xv() const {
  if (leading_index() < 0 || detected_photon_index() < 0 ||
      target_index() < 0) {
    return 0.;
  }
  return (detected_Q2() + leading().mass2()) /
         (2 * target().mass() * detected_nu());
}
inline double lp_gamma_event::detected_Q2plusM2() const {
  if (detected_leading_index() < 0 || detected_photon_index() < 0) {
    return 0.;
  }
  return detected_Q2() + detected_leading().mass2();
}

} // namespace pcsim

#endif
