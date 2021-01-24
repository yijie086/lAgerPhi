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

#include "jpacPhoto_pomeron.hh"
#include <TMath.h>
#include <lager/core/logger.hh>
#include <lager/gen/initial/target_gen.hh>
#include <lager/physics/decay.hh>
#include <lager/physics/kinematics.hh>
#include <lager/physics/photon.hh>
#include <lager/physics/vm.hh>

namespace lager {
namespace lA {

// =============================================================================
// Constructor for lA::jpacPhoto_pomeron
// =============================================================================
jpacPhoto_pomeron::jpacPhoto_pomeron(const configuration& cf,
                                     const string_path& path,
                                     std::shared_ptr<TRandom> r)
    : base_type{r}
    , recoil_{cf.get<std::string>(path / "recoil_type")}
    , vm_{cf.get<std::string>(path / "vm_type")}
    , decay_lplus_{static_cast<pdg_id>(
          -abs(cf.get<int>(path / "vm_decay_lepton_type", 11)))}
    , decay_lminus_{static_cast<pdg_id>(
          abs(cf.get<int>(path / "vm_decay_lepton_type", 11)))}
    , regge_inter_{cf.get<double>(path / "regge_inter")}
    , regge_slope_{cf.get<double>(path / "regge_slope")}
    , photo_norm_{cf.get<double>(path / "photo_norm")}
    , photo_slope_{cf.get<double>(path / "photo_slope")}
    , dipole_n_{cf.get<double>(path / "dipole_n")}
    , reaction_{init_reaction()}
    , regge_{1, regge_inter_, regge_slope_, "pomeron"}
    , ampl_{init_ampl()}
    , max_t_range_{calc_max_t_range(cf)}
    , max_{calc_max_xsec(cf)} {
  LOG_INFO("jpacPhoto_pomeron", "t range [GeV^2]: [" +
                                    std::to_string(max_t_range_.min) + ", " +
                                    std::to_string(max_t_range_.max) + "]");
  LOG_INFO("jpacPhoto_pomeron",
           "regge trajectory intercept: " + std::to_string(regge_inter_));
  LOG_INFO("jpacPhoto_pomeron",
           "regge trajectory slope: " + std::to_string(regge_slope_));
  LOG_INFO("jpacPhoto_pomeron",
           "photo norm parameter [1/GeV^2]: " + std::to_string(photo_norm_));
  LOG_INFO("jpacPhoto_pomeron",
           "photo slope parameter [1/GeV^2]: " + std::to_string(photo_slope_));
  LOG_INFO("jpacPhoto_pomeron",
           "'Dipole' FF power: " + std::to_string(dipole_n_));
  LOG_INFO("jpacPhoto_pomeron", "VM: " + std::string(vm_.pdg()->GetName()));
  LOG_INFO("jpacPhoto_pomeron",
           "recoil: " + std::string(recoil_.pdg()->GetName()));
}
std::unique_ptr<jpacPhoto::reaction_kinematics>
jpacPhoto_pomeron::init_reaction() const {
  auto reaction = std::make_unique<jpacPhoto::reaction_kinematics>(vm_.mass());
  reaction->set_JP(1, -1);
  return reaction;
}
jpacPhoto::pomeron_exchange jpacPhoto_pomeron::init_ampl() {
  jpacPhoto::pomeron_exchange ampl{reaction_.get(), &regge_, false,
                                   "t-channel"};
  ampl.set_params({photo_norm_, photo_slope_});
  return ampl;
}

lA_event jpacPhoto_pomeron::generate(const lA_data& initial) {
  // generate a mass() in case of non-zero width, initialize the particles
  particle vm = {vm_.type(), rng()};
  particle recoil = {recoil_.type(), rng()};

  // shortcuts
  const auto& gamma = initial.photon();
  const auto& target = initial.target();

  // check if enough energy available
  if (gamma.W2() < threshold2(vm, recoil)) {
    LOG_JUNK("jpacPhoto_pomeron", "Not enough phase space available - W2: " +
                                      std::to_string(gamma.W2()) + " < " +
                                      std::to_string(threshold2(vm, recoil)));
    return lA_event{0.};
  }

  // generate a phase space point
  const double t = rng()->Uniform(max_t_range_.min, max_t_range_.max);

  LOG_JUNK("jpacPhoto_pomeron", "t: " + std::to_string(t));

  // check if kinematically allowed
  if (physics::t_range(gamma.W2(), gamma.Q2(), target.particle().mass(),
                       vm.mass(), recoil.mass())
          .excludes(t)) {
    LOG_JUNK("jpacPhoto_pomeron", "t outside of the allowed range for this W2")
    return lA_event{0.};
  }

  // evaluate the cross section
  const double xs_dipole = dipole(gamma.Q2());
  const double xs_photo = dsigma_dt(gamma.W2(), t);
  const double xs = xs_dipole * xs_photo;

  LOG_JUNK("jpacPhoto_pomeron",
           "xsec: " + std::to_string(xs_photo) + " < " + std::to_string(max_));
  LOG_JUNK("jpacPhoto_pomeron", "dipole: " + std::to_string(xs_dipole));

  // return a new VM event
  return make_event(initial, t, vm, recoil, xs);
}

// =============================================================================
// jpacPhoto_pomeron::calc_max_xsec(cf)
//
// Utility function for the generator initialization
//
// max cross section as defined by the beam/target and phase-space settings
//
// The cross section is maximum when
//  * W2 is maximum.
//  * no t-requirement (we throw directly in exp(bt)
//  * no Q2 requirement:
//    - irrelevant for photo-production
//    - for lepto-production, the (1 + epsilon * R) * dipole factor does not
//      change the value for the fully differential cross section maximum at
//      Q2min
//  * note that we can be certain about these statements, as the program will
//    exit with an error if the cross section maximum were ever violated
// =============================================================================
double jpacPhoto_pomeron::calc_max_xsec(const configuration& cf) /*const */ {
  // get the extreme beam parameters (where the photon carries all of the
  // lepton beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                        cf.get<double>("beam/lepton/energy")};
  const particle target{initial::estimated_target(cf)};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W
  const double W2max = opt_W_range ? fmin(opt_W_range->max * opt_W_range->max,
                                          (photon.p() + target.p()).M2())
                                   : (photon.p() + target.p()).M2();
  return dsigma_dt(W2max, max_t_range_.max) * 1.0001;
} // namespace lA

// =============================================================================
// jpacPhoto_pomeron::calc_max_t_range(cf)
//
// Utility function for the generator initialization
//
// max t-range occurs for:
//  * photon that carries all of the beam energy (or when we have reached the
//    user-defined maximum value of W max)
//  * maximum Q2 for the given W (or zero for real photons) for the lower t
//  bound
//  * Q2 = 0 for the upper t bound (tmin)
// =============================================================================
interval<double>
jpacPhoto_pomeron::calc_max_t_range(const configuration& cf) const {
  // hack
  // get the extreme beam parameters (where the photon carries all of the
  // lepton beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                        cf.get<double>("beam/lepton/energy")};
  const particle target{initial::estimated_target(cf)};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W (where the t-range is the largest)
  const double W2max = opt_W_range ? fmin(opt_W_range->max * opt_W_range->max,
                                          (photon.p() + target.p()).M2())
                                   : (photon.p() + target.p()).M2();
  // some shortcuts
  const double M_nu = (photon.p()).Dot(target.p());
  const double Q2max = target.mass2() + 2 * M_nu - W2max;
  const double Q2min = 0.;

  // calculate the corresponding t range
  // in case of particles with non-zero width, we use M - 4 x sigma
  //
  // tlim1: the correct lower t limit (tmax)
  const auto tlim1 =
      physics::t_range(W2max, (Q2max > 1e-10 ? Q2max : 0.), target.mass(),
                       vm_.pole_mass() - vm_.width() * 4.,
                       recoil_.pole_mass() - recoil_.width() * 4);
  // tlim2: the correct upper t limit (tmin)
  const auto tlim2 = physics::t_range(
      W2max, 0, target.mass(), vm_.pole_mass() - vm_.width() * 4.,
      recoil_.pole_mass() - recoil_.width() * 4);
  const auto tlim = interval<double>(tlim1.min, tlim2.max);
  return tlim;
}
// =============================================================================
// jpacPhoto_pomeron::dsigma_dt()
// jpacPhoto_pomeron::dipole()
//
// Utility functions to calculate the cross section components
// =============================================================================
double jpacPhoto_pomeron::dsigma_dt(const double s, const double t) /*const */ {
  return ampl_.differential_xsection(s, t);
}
double jpacPhoto_pomeron::dipole(const double Q2) const {
  return physics::dipole_ff_vm(Q2, vm_.mass(), dipole_n_);
}
// =============================================================================
// jpacPhoto_pomeron::threshold2()
//
// utility function returns the production threshold squared for the chosen
// particles. Important to re-calculate in case of a particle with non-zero
// widht.
// =============================================================================
double jpacPhoto_pomeron::threshold2(const particle& vm,
                                     const particle& recoil) const {
  return recoil.mass2() + vm.mass2() + 2 * vm.mass() * recoil.mass();
}

// =============================================================================
// create the lA_event data, calculates the final state four-vectors in
// the lab-frame
// =============================================================================
lA_event jpacPhoto_pomeron::make_event(const lA_data& initial, const double t,
                                       particle vm, particle X,
                                       const double xs) {
  const auto& gamma = initial.photon();
  const auto& target = initial.target();

  lA_event e{initial, xs, 1., 0 /* R */};

  // utility shortcuts
  const double W2 = gamma.W2();
  const double W = sqrt(W2);
  const double Q2 = gamma.Q2();
  const double y = gamma.y();
  const double nu = gamma.nu();
  const double Mt2 = target.particle().mass2();
  const double Mr2 = X.mass2();
  const double Mv2 = vm.mass2();

  // create our final state particles in the CM frame
  // energies and momenta
  const double Et_cm = (W2 + Q2 + Mt2) / (2. * W);
  const double Pt_cm = sqrt(Et_cm * Et_cm - Mt2);
  const double Er_cm = (W2 - Mv2 + Mr2) / (2. * W);
  const double Pr_cm = sqrt(Er_cm * Er_cm - Mr2);
  const double Ev_cm = (W2 + Mv2 - Mr2) / (2. * W);
  const double Pv_cm = sqrt(Ev_cm * Ev_cm - Mv2);

  // get the VM and recoil CM 4-vectors
  // NOTE:
  //  * theta is the change in angle from the initial state, don't forget the
  //    target originally was flying backwards (already has theta=pi)
  // calculate the scattering angle theta
  const double ctheta_cm =
      (t + 2 * Et_cm * Er_cm - Mt2 - Mr2) / (2 * Pt_cm * Pr_cm);
  const double theta_cm = std::acos(ctheta_cm);
  const double phi_cm = rng()->Uniform(0, TMath::TwoPi());

  // create the CM 4-vectors
  const particle::Polar3DVector p3_v{Pv_cm, theta_cm, phi_cm};
  vm.p() = {p3_v.X(), p3_v.Y(), p3_v.Z(), Ev_cm};
  const particle::Polar3DVector p3_r{Pr_cm, theta_cm + TMath::Pi(), phi_cm};
  X.p() = {p3_r.X(), p3_r.Y(), p3_r.Z(), Er_cm};

  // VM decay particles
  std::pair<particle, particle> decay_products{{decay_lplus_.type()},
                                               {decay_lminus_.type()}};
  auto& d0 = decay_products.first;
  auto& d1 = decay_products.second;
  {
    // Also create CM 4-vector in the production plane (phi==0), needed to
    // properly handle the VM decay. Note that this is just a rotated version
    // of the VM 4-vector
    particle vm_pplane = vm;
    const particle::Polar3DVector p3_v_pplane{Pv_cm, theta_cm, 0};
    vm_pplane.p() = {p3_v_pplane.X(), p3_v_pplane.Y(), p3_v_pplane.Z(), Ev_cm};
    // deprecated CM info of decay products
    std::pair<particle, particle> decay_products_cm{
        {decay_lplus_.type(), particle::status_code::INFO_PARENT_CM},
        {decay_lminus_.type(), particle::status_code::INFO_PARENT_CM}};
    // get the real part of the SDME for this event
    const double r00 = std::real(ampl_.SDME(0, 0, 0, W2, t));
    const double r10 = std::real(ampl_.SDME(0, 1, 0, W2, t));
    const double r1m1 = std::real(ampl_.SDME(0, 1, -1, W2, t));
    LOG_JUNK2("jpacPhoto_pomeron", "SDME Values: r00: " + std::to_string(r00) +
                                       ", r10: " + std::to_string(r10) +
                                       ", r1-1: " + std::to_string(r1m1));
    auto decay_distribution = [&](double cth, double phi) {
      return physics::vm_decay_fermions(cth, phi, r00, r10, r1m1);
    };
    const auto [ctheta, phi] =
        rand_f({-1, 1}, {0., TMath::TwoPi()}, decay_distribution, 2.0);
    const double theta = acos(ctheta);
    physics::decay_2body(vm_pplane, theta, phi, decay_products,
                         decay_products_cm);
    // Rotate back around the virtual photon direction to the true CM frame
    const particle::Polar3DVector v0{d0.momentum(), d0.theta(),
                                     d0.phi() + phi_cm};
    const particle::Polar3DVector v1{d1.momentum(), d1.theta(),
                                     d1.phi() + phi_cm};
    d0.p() = {v0.X(), v0.Y(), v0.Z(), d0.energy()};
    d1.p() = {v1.X(), v1.Y(), v1.Z(), d1.energy()};
  }

  // calculate some necessary boost and rotation vectors
  particle targ = target.particle();
  particle phot = gamma.particle();
  // lab frame to target rest frame (trf)
  particle::Boost boost_to_trf{targ.p().BoostToCM()};
  targ.boost(boost_to_trf);
  phot.boost(boost_to_trf);
  // describe photon in a target rest frame where the photon moves along the
  // z-axis
  particle::XYZTVector p_phot_z{0, 0, phot.momentum(), phot.energy()};
  // boost to CM frame from this rotated trf
  particle::Boost boost_from_cm{-((targ.p() + p_phot_z).BoostToCM())};

  // go to lab frame
  // 1. CM -> rotated TRF
  vm.boost(boost_from_cm);
  X.boost(boost_from_cm);
  d0.boost(boost_from_cm);
  d1.boost(boost_from_cm);
  // 2. rotated TRF -> TRF
  vm.rotate_uz(phot.p());
  X.rotate_uz(phot.p());
  d0.rotate_uz(phot.p());
  d1.rotate_uz(phot.p());
  // 3. TRF -> lab
  vm.boost(boost_to_trf.Inverse());
  X.boost(boost_to_trf.Inverse());
  d0.boost(boost_to_trf.Inverse());
  d1.boost(boost_to_trf.Inverse());

  // update VM particle status
  vm.update_status(particle::status_code::UNSTABLE_RADCOR_ONLY);

  // set vertex info
  vm.vertex() = phot.vertex();
  X.vertex() = phot.vertex();
  d0.vertex() = phot.vertex();
  d1.vertex() = phot.vertex();

  // add to the event
  int vm_idx = e.add_leading(vm, e.photon_index(), e.target_index());
  e.add_recoil(X, e.photon_index(), e.target_index());
  e.add_daughter(decay_products, vm_idx);

  // all done!
  return e;
}

} // namespace lA
} // namespace lager
