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

#include "oleksii_jpsi_bh.hh"
#include "oleksii_bh_impl.hh"
#include "oleksii_total_impl.hh"
#include <TF3.h>
#include <TMath.h>
#include <lager/core/assert.hh>
#include <lager/core/logger.hh>
#include <lager/physics/decay.hh>
#include <lager/physics/kinematics.hh>
#include <lager/physics/photon.hh>
#include <lager/physics/vm.hh>

namespace lager {
namespace lA {

// =============================================================================
// Constructor for lA::oleksii_jpsi_bh
// =============================================================================
oleksii_jpsi_bh::oleksii_jpsi_bh(const configuration& cf,
                                 const string_path& path,
                                 std::shared_ptr<TRandom> r)
    : base_type{r}
    , recoil_{pdg_id::p}
    , vm_{pdg_id::J_psi}
    , max_t_range_{calc_max_t_range(cf)}
    , Mll_range_{cf.get_range<double>(path / "Mll_range")}
    , max_{cf.get<double>(path / "max_cross_section")}
    , T_0_{cf.get<double>(path / "T_0")}
    , theta_range_{cf.get_range<double>(path / "theta_range") *
                   TMath::DegToRad()}
    , p_range_{cf.get_range<double>(path / "p_range")} {
  LOG_INFO("oleksii_jpsi_bh", "t range [GeV^2]: [" +
                                  std::to_string(max_t_range_.min) + ", " +
                                  std::to_string(max_t_range_.max) + "]");
  LOG_INFO("oleksii_jpsi_bh", "Mll range [GeV]: [" +
                                  std::to_string(Mll_range_.min) + ", " +
                                  std::to_string(Mll_range_.max) + "]");
  LOG_INFO("oleksii_jpsi_bh",
           "Maximum cross section set to: " + std::to_string(max_));
  LOG_INFO("oleksii_jpsi_bh",
           "Subtraction constant T_0: " + std::to_string(T_0_));
  LOG_INFO("oleksii_jpsi_bh",
           "Theta acceptance [deg.]: [" +
               std::to_string(theta_range_.min * TMath::RadToDeg()) + ", " +
               std::to_string(theta_range_.max * TMath::RadToDeg()) + "]");
  LOG_INFO("oleksii_jpsi_bh", "P acceptance [GeV]: [" +
                                  std::to_string(p_range_.min) + ", " +
                                  std::to_string(p_range_.max) + "]");
}

lA_event oleksii_jpsi_bh::generate(const lA_data& initial) {

  // particle vm = {};
  particle recoil = {pdg_id::p};

  // shortcuts
  const auto& gamma = initial.photon();
  const auto& target = initial.target();

  // throw error if we are running electro-production
  tassert(gamma.Q2() < 1,
          "oleksii_jpsi_bh is a photo-production only generator!");

  // generate a Mll point
  const double Mll2 = rng()->Uniform(Mll_range_.min * Mll_range_.min,
                                     Mll_range_.max * Mll_range_.max);
  const double Mll = sqrt(Mll2);

  LOG_JUNK("oleksii_jpsi_bh", "Mll2: " + std::to_string(Mll2));
  LOG_JUNK("oleksii_jpsi_bh", "(Mll: " + std::to_string(Mll) + ")");

  // create our "vm"
  particle vm = {pdg_id::J_psi, {0, 0, 0, Mll}};

  // check if enough energy available
  if (gamma.W2() < threshold2(vm, recoil)) {
    LOG_JUNK("oleksii_jpsi_bh", "Not enough phase space available - W2: " +
                                    std::to_string(gamma.W2()) + " < " +
                                    std::to_string(threshold2(vm, recoil)));
    return lA_event{0.};
  }

  const double t = std::log(rng()->Uniform(std::exp(1.13 * max_t_range_.min),
                                           std::exp(1.13 * max_t_range_.max))) /
                   1.13;
  const double jacobian = exp(-1.13 * t) / 1.13;
  // further generate this phase space point
  // const double t = std::log(rng()->Uniform(std::exp(5 * max_t_range_.min),
  //                                         std::exp(5 * max_t_range_.max)))
  //                                         /
  //                 5.;
  // const double t = rng()->Uniform(max_t_range_.min, max_t_range_.max);
  // const double thetaCM = acos(rng()->Uniform(-0.766, -0.707));
  const double thetaCM = acos(rng()->Uniform(-1., 1.));
  const double phiCM = rng()->Uniform(0, TMath::TwoPi());

  LOG_JUNK("oleksii_jpsi_bh", "t: " + std::to_string(t));
  LOG_JUNK("oleksii_jpsi_bh", "thetaCM: " + std::to_string(thetaCM));
  LOG_JUNK("oleksii_jpsi_bh", "phiCM: " + std::to_string(phiCM));

  // check if kinematically allowed
  auto tlim = physics::t_range(gamma.W2(), gamma.Q2(), target.particle().mass(),
                               vm.mass(), recoil.mass());
  if (tlim.excludes(t)) {
    LOG_JUNK("oleksii_jpsi_bh", "t outside of the allowed range for this W2");
    return lA_event{0.};
  }

  // invariant definition of Egamma
  const double Egamma =
      (gamma.particle().p()).Dot(target.particle().p()) / target.particle().mass();

  // evaluate the cross section
  const double xs =
      oleksii_total_impl::calc_xsec(t, Mll2, Egamma, thetaCM, phiCM, T_0_);

  LOG_JUNK("oleksii_jpsi_bh",
           "xsec: " + std::to_string(xs) + " < " + std::to_string(max_));

  // return a new VM event
  return make_event(initial, t, vm, recoil, xs * jacobian, thetaCM, phiCM);
}

// =============================================================================
// oleksii_jpsi_bh::calc_max_t_range(cf)
//
// Utility function for the generator initialization
//
// max t-range occurs for:
//  * photon that carries all of the beam energy (or when we have reached the
//    user-defined maximum value of W max)
//  * maximum Q2 for the given W (or zero for real photons)
// =============================================================================
interval<double>
oleksii_jpsi_bh::calc_max_t_range(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the
  // lepton beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                        cf.get<double>("beam/lepton/energy")};
  const particle target{
      cf.get<std::string>("beam/ion/particle_type"),
      cf.get_vector3<particle::XYZVector>("beam/ion/dir"),
      cf.get<double>("beam/ion/energy")};
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
  // const auto tlim = interval<double>(tlim1.min, tlim2.max*1.01);
  // SJJ hardcoded tmax
  const auto tlim = interval<double>(tlim2.max - 4.0, tlim2.max - .0001);
  return tlim;
  // HARD CODED
  // return {-1.0, -0.95};
}

// =============================================================================
// oleksii_jpsi_bh::threshold2()
//
// utility function returns the production threshold squared for the chosen
// particles. Important to re-calculate in case of a particle with non-zero
// widht.
// =============================================================================
double oleksii_jpsi_bh::threshold2(const particle& vm,
                                   const particle& recoil) const {

  return recoil.mass2() + vm.mass2() + 2 * vm.mass() * recoil.mass();
}

// =============================================================================
// create the lA_event dataf, calculates the final state four-vectors in
// the lab-frame
// =============================================================================
lA_event oleksii_jpsi_bh::make_event(const lA_data& initial,
                                           const double t, particle vm,
                                           particle X, const double xs,
                                           const double thetaCM_el,
                                           const double phiCM_el) const {
  const auto& gamma = initial.photon();
  const auto& target = initial.target();
  lA_event e{initial, xs, 1., 0};

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

  // first work in the CM frame to handle the phi dependence correctly

  // 'decay' the 'VM' into our lepton pair
  std::pair<particle, particle> decay_products{{pdg_id::e_minus},
                                               {pdg_id::e_plus}};
  std::pair<particle, particle> decay_products_cm{
      {pdg_id::e_minus, particle::status_code::INFO_PARENT_CM},
      {pdg_id::e_plus, particle::status_code::INFO_PARENT_CM}};

  physics::decay_2body(vm, thetaCM_el, phiCM_el, decay_products,
                       decay_products_cm);

  // calculate some necessary boost and rotation vectors to go to the lab
  // frame
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
  decay_products.first.boost(boost_from_cm);
  decay_products.second.boost(boost_from_cm);
  // 2. rotated TRF -> TRF
  vm.rotate_uz(phot.p());
  X.rotate_uz(phot.p());
  decay_products.first.rotate_uz(phot.p());
  decay_products.second.rotate_uz(phot.p());
  // 3. TRF -> lab
  vm.boost(boost_to_trf.Inverse());
  X.boost(boost_to_trf.Inverse());
  decay_products.first.boost(boost_to_trf.Inverse());
  decay_products.second.boost(boost_to_trf.Inverse());

  // update VM particle status
  vm.update_status(particle::status_code::DECAYED);

  // add vertex
  vm.vertex() = phot.vertex();
  X.vertex() = phot.vertex();
  decay_products.first.vertex() = phot.vertex();
  decay_products.second.vertex() = phot.vertex();

  // add to the event
  int vm_idx = e.add_leading(vm, e.photon_index(), e.target_index());
  e.add_recoil(X, e.photon_index(), e.target_index());

  e.add_daughter(decay_products, vm_idx);
  decay_products_cm.first.add_parent(vm_idx);
  decay_products_cm.second.add_parent(vm_idx);
  e.add_particle(decay_products_cm);

  // momentum and polar angle cut
  const double pos_mom = e[7].momentum();
  const double el_mom = e[6].momentum();
  const double pos_theta = e[7].p().Theta();
  const double el_theta = e[6].p().Theta();

  if (p_range_.excludes(pos_mom) || p_range_.excludes(el_mom) ||
      theta_range_.excludes(pos_theta) || theta_range_.excludes(el_theta)) {
    return lA_event{0.};
  } else {
    ; // do nothing
  }

  // all done!
  return e;
}

} // namespace lA
} // namespace lager
