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

#ifndef LAGER_GEN_LA_EVENT_LOADED
#define LAGER_GEN_LA_EVENT_LOADED

#include <fstream>
#include <lager/core/event.hh>
#include <lager/core/generator.hh>
#include <lager/gen/initial/data.hh>

namespace lager {

// =============================================================================
// Initial reaction data for a lepton-proton reaction mediated by a virtual
// photon.
// =============================================================================
class lA_data : public generator_data {
public:
  lA_data(const double xs = 0) : generator_data{xs} {}
  lA_data(const initial::beam& l, const initial::beam& p,
          const initial::target& targ, const initial::photon& gamma);
  lA_data(const lA_data&) = default;
  lA_data& operator=(const lA_data&) = default;

  initial::beam& lepton() { return lepton_; }
  const initial::beam& lepton() const { return lepton_; }
  initial::beam& ion() { return ion_; }
  const initial::beam& ion() const { return ion_; }
  initial::target& target() { return target_; }
  const initial::target& target() const { return target_; }
  initial::photon& photon() { return photon_; }
  const initial::photon& photon() const { return photon_; }

private:
  initial::beam lepton_;
  initial::beam ion_;
  initial::target target_;
  initial::photon photon_;
};

// =============================================================================
// full lA event
// =============================================================================
class lA_event : public event {
public:
  lA_event(const lA_event&) = default;
  lA_event& operator=(const lA_event&) = default;
  explicit lA_event(const double xs = 1., const double w = 1.,
                    const double R = 0., const double epsilon = 0.);
  lA_event(const lA_data& initial, const double xs, const double w = 1.,
           const double R = 0.);

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
  // add target and possible spectator particles
  int add_target(const initial::target& target);
  // add photon and scat
  int add_photon(const initial::photon& gamma);

  // add leading particle
  int add_leading(const particle& part, const int parent1,
                  const int parent2 = -1);
  int add_leading(const particle& part);
  // add recoil particle
  int add_recoil(const particle& part, const int parent1,
                 const int parent2 = -1);
  int add_recoil(const particle& part);

  // get photon, scat, leading and recoil info
  particle& target();
  const particle& target() const;
  particle& photon();
  const particle& photon() const;
  particle& scat();
  const particle& scat() const;
  particle& leading();
  const particle& leading() const;
  particle& recoil();
  const particle& recoil() const;

  int target_index() const { return target_index_; }
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
  int target_index_{-1};
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
// lA_out custom event output
// =============================================================================
class lA_out : public event_out {
public:
  lA_out(std::shared_ptr<TFile> f, std::unique_ptr<std::ofstream> olund,
         std::unique_ptr<std::ofstream> osimc, const std::string& name);

  void push(const lA_event& e);
  void push(const std::vector<lA_event>& e);

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
  int16_t target_index_;
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
// LA_DATA IMPLEMENTATION
// =============================================================================
inline lA_data::lA_data(const initial::beam& l, const initial::beam& A,
                        const initial::target& t, const initial::photon& gamma)
    : generator_data{gamma.cross_section()}
    , lepton_{l}
    , ion_{A}
    , target_{t}
    , photon_{gamma} {}

// =============================================================================
// LA_EVENT IMPLEMENTATION
// =============================================================================

// constructors
inline lA_event::lA_event(const double xs, const double w, const double R,
                          const double epsilon)
    : event{xs, w}, R_{R}, epsilon_{epsilon} {}
inline lA_event::lA_event(const lA_data& initial, const double xs,
                          const double w, const double R)
    : event{xs, w}, R_{R}, epsilon_{initial.photon().epsilon()} {
  add_ibeam(initial.lepton().particle());
  add_tbeam(initial.ion().particle());
  // thes step also applies the target/photon cross section!
  add_target(initial.target());
  add_photon(initial.photon());
}

// set particles
inline int lA_event::add_target(const initial::target& targ) {
  target_index_ = add_daughter(targ.particle(), tbeam_index());
  for (const auto& remn : targ.remnant()) {
    add_daughter(remn, tbeam_index());
  }

  update_cross_section(targ.cross_section());

  return target_index();
}
inline int lA_event::add_photon(const initial::photon& gamma) {
  photon_index_ = add_daughter(gamma.particle(), ibeam_index());
  scat_index_ = add_daughter(gamma.scat(), ibeam_index());

  W_ = sqrt(gamma.W2());
  Q2_ = gamma.Q2();
  nu_ = gamma.nu();
  x_ = gamma.x();
  y_ = gamma.y();

  update_cross_section(gamma.cross_section());

  return photon_index();
}
inline int lA_event::add_leading(const particle& part, const int parent1,
                                 const int parent2) {
  leading_index_ = add_daughter(part, parent1, parent2);
  return leading_index();
}
inline int lA_event::add_leading(const particle& part) {
  return add_leading(part, photon_index(), target_index());
}
inline int lA_event::add_recoil(const particle& part, const int parent1,
                                const int parent2) {
  recoil_index_ = add_daughter(part, parent1, parent2);
  return recoil_index();
}
inline int lA_event::add_recoil(const particle& part) {
  return add_recoil(part, photon_index(), target_index());
}

// get target, photon, scat, leading and recoil info
inline particle& lA_event::target() {
  tassert(target_index_ >= 0, "trying to access target data, but no target "
                              "data present in the event");
  return part(target_index());
}
inline const particle& lA_event::target() const {
  tassert(target_index_ >= 0, "trying to access target data, but no target "
                              "data present in the event");
  return part(target_index());
}
inline particle& lA_event::photon() {
  tassert(photon_index_ >= 0, "trying to access photon data, but no photon "
                              "data present in the event");
  return part(photon_index());
}
inline const particle& lA_event::photon() const {
  tassert(photon_index_ >= 0, "trying to access photon data, but no photon "
                              "data present in the event");
  return part(photon_index());
}
inline particle& lA_event::scat() {
  tassert(scat_index_ >= 0, "trying to access scat data, but no scat "
                            "data present in the event");
  return part(scat_index());
}
inline const particle& lA_event::scat() const {
  tassert(scat_index_ >= 0, "trying to access scat data, but no scat "
                            "data present in the event");
  return part(scat_index());
}
inline particle& lA_event::leading() {
  tassert(leading_index_ >= 0, "trying to access leading data, but no leading "
                               "data present in the event");
  return part(leading_index());
}
inline const particle& lA_event::leading() const {
  tassert(leading_index_ >= 0, "trying to access leading data, but no leading "
                               "data present in the event");
  return part(leading_index());
}
inline particle& lA_event::recoil() {
  tassert(recoil_index_ >= 0, "trying to access recoil data, but no recoil "
                              "data present in the event");
  return part(recoil_index());
}
inline const particle& lA_event::recoil() const {
  tassert(recoil_index_ >= 0, "trying to access recoil data, but no recoil "
                              "data present in the event");
  return part(recoil_index());
}
// detected particle, also set the relevant detected index
inline void lA_event::update_detected_index(const int index) {
  const auto& dp = detected(index);
  if (dp.generated().index() == scat_index()) {
    LOG_JUNK2("lA_event", "Found reconstructed scattered lepton");
    rc_scat_index_ = index;
    // also reconstruct the photon
    // status set to -1 for purely reconstucted particle
    if (ibeam_index() >= 0 && photon_index() >= 0) {
      LOG_JUNK2("lA_event", "Adding reconstructed photon data");
      rc_photon_index_ =
          add_detected({photon(), particle().p() - detected_scat().p(), -1});
    } else {
      LOG_JUNK2("lA_event", "Unable to add reconstructed photon data");
    }
  } else if (dp.generated().index() == photon_index()) {
    LOG_JUNK2("lA_event", "Found reconstructed scattered photon");
    rc_photon_index_ = index;
  } else if (dp.generated().index() == leading_index()) {
    LOG_JUNK2("lA_event", "Found reconstructed leading particle");
    rc_leading_index_ = index;
  } else if (dp.generated().index() == recoil_index()) {
    LOG_JUNK2("lA_event", "Found reconstructed recoil particle");
    rc_recoil_index_ = index;
  }
}

inline detected_particle& lA_event::detected_photon() {
  tassert(detected_photon_index() >= 0,
          "trying to access detected_photon data, but no detected_photon "
          "data present in the event");
  return detected(detected_photon_index());
}
inline const detected_particle& lA_event::detected_photon() const {
  tassert(detected_photon_index() >= 0,
          "trying to access detected_photon data, but no detected_photon "
          "data present in the event");
  return detected(detected_photon_index());
}
inline detected_particle& lA_event::detected_scat() {
  tassert(detected_scat_index() >= 0,
          "trying to access detected_scat data, but no detected_scat "
          "data present in the event");
  return detected(detected_scat_index());
}
inline const detected_particle& lA_event::detected_scat() const {
  tassert(detected_scat_index() >= 0,
          "trying to access detected_scat data, but no detected_scat "
          "data present in the event");
  return detected(detected_scat_index());
}
inline detected_particle& lA_event::detected_leading() {
  tassert(detected_leading_index() >= 0,
          "trying to access detected_leading data, but no detected_leading "
          "data present in the event");
  return detected(detected_leading_index());
}
inline const detected_particle& lA_event::detected_leading() const {
  tassert(detected_leading_index() >= 0,
          "trying to access detected_leading data, but no detected_leading "
          "data present in the event");
  return detected(detected_leading_index());
}
inline detected_particle& lA_event::detected_recoil() {
  tassert(detected_recoil_index() >= 0,
          "trying to access detected_recoil data, but no detected_recoil "
          "data present in the event");
  return detected(detected_recoil_index());
}
inline const detected_particle& lA_event::detected_recoil() const {
  tassert(detected_recoil_index() >= 0,
          "trying to access detected_recoil data, but no detected_recoil "
          "data present in the event");
  return detected(detected_recoil_index());
}

// leading produced particle kinematics
inline double lA_event::t() const {
  if (leading_index() >= 0 && photon_index() >= 0) {
    return (leading().p() - photon().p()).M2();
  } else if (target_index() >= 0 && recoil_index() >= 0) {
    return (target().p() - recoil().p()).M2();
  }
  return 0.;
}
inline double lA_event::xv() const {
  if (leading_index() < 0 || photon_index() < 0 || target_index() < 0) {
    return 0.;
  }
  return (Q2() + leading().mass2()) / (2 * target().mass() * nu());
}
inline double lA_event::Q2plusM2() const {
  if (leading_index() < 0 || photon_index() < 0) {
    return 0.;
  }
  return Q2() + leading().mass2();
}

// detected kinematics
inline double lA_event::detected_W() const { return sqrt(detected_W2()); }
inline double lA_event::detected_W2() const {
  if (detected_photon_index() < 0 || tbeam_index() < 0) {
    return 0.;
  }
  return (detected_photon().p() + tbeam().p()).M2();
}
inline double lA_event::detected_Q2() const {
  if (detected_photon_index() < 0) {
    return 0.;
  }
  return -detected_photon().p().M2();
}
inline double lA_event::detected_nu() const {
  if (detected_photon_index() < 0 || tbeam_index() < 0) {
    return 0.;
  }
  return (detected_photon().p()).Dot(tbeam().p()) / target().mass();
}
inline double lA_event::detected_x() const {
  if (detected_photon_index() < 0 || tbeam_index() < 0) {
    return 0.;
  }
  return detected_Q2() / (2 * tbeam().mass() * detected_nu());
}
inline double lA_event::detected_y() const {
  if (detected_photon_index() < 0 || tbeam_index() < 0 || ibeam_index() < 0) {
    return 0.;
  }
  return (detected_photon().p()).Dot(tbeam().p()) /
         (particle().p()).Dot(tbeam().p());
}
inline double lA_event::detected_t() const {
  if (detected_photon_index() < 0 || detected_leading_index() < 0) {
    return 0.;
  }
  return (detected_leading().p() - detected_photon().p()).M2();
}
inline double lA_event::detected_xv() const {
  if (leading_index() < 0 || detected_photon_index() < 0 || tbeam_index() < 0) {
    return 0.;
  }
  return (detected_Q2() + leading().mass2()) /
         (2 * tbeam().mass() * detected_nu());
}
inline double lA_event::detected_Q2plusM2() const {
  if (detected_leading_index() < 0 || detected_photon_index() < 0) {
    return 0.;
  }
  return detected_Q2() + detected_leading().mass2();
}

} // namespace lager

#endif
