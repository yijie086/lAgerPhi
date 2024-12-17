// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
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

#include "phi_hatta.hh"
#include <TMath.h>
//#include <lager/core/factory.hh>
#include <lager/core/logger.hh>
#include <lager/gen/initial/target_gen.hh>
#include <lager/physics/kinematics.hh>
#include <lager/physics/photon.hh>
#include <lager/physics/vm.hh>

namespace lager {
namespace lA {

// =============================================================================
// Constructor for lA::phi_hatta
// =============================================================================
phi_hatta::phi_hatta(const configuration& cf, const string_path& path,
                       std::shared_ptr<TRandom> r)
    : base_type{r}
    , recoil_{cf.get<std::string>(path / "recoil_type")}
    , vm_{cf.get<std::string>(path / "vm_type")}
    , alpha_1_{cf.get<double>(path / "alpha_1")}
    , alpha_2_{cf.get<double>(path / "alpha_2")}
    , alpha_3_{cf.get<double>(path / "alpha_3")}
    , nu_T_{cf.get<double>(path / "nu_T")}
    , c_R_{cf.get<double>(path / "c_R")}
    , max_t_range_{calc_max_t_range(cf, path)}
    , max_{calc_max_xsec(cf)} {
  LOG_INFO("phi_hatta", "t range [GeV^2]: [" +
                             std::to_string(max_t_range_.min) + ", " +
                             std::to_string(max_t_range_.max) + "]");
  LOG_INFO("phi_hatta", "alpha_1 parameter: " + std::to_string(alpha_1_));
  LOG_INFO("phi_hatta", "alpha_2 parameter: " + std::to_string(alpha_2_));
  LOG_INFO("phi_hatta", "alpha_3 parameter: " + std::to_string(alpha_3_));
  LOG_INFO("phi_hatta", "nu_T parameter: " + std::to_string(nu_T_));
  LOG_INFO("phi_hatta", "c_R parameter: " + std::to_string(c_R_));
  const std::string ff = cf.get<std::string>(path / "ff" / "function");
  LOG_INFO("phi_hatta", "FF function: " + ff);
  if (ff == "exp") {
    const double B0 = cf.get<double>(path / "ff" / "B0");
    const double alphaP = cf.get<double>(path / "ff" / "alphaP");
    LOG_INFO("phi_hatta", "FF B0 parameter" + std::to_string(B0));
    LOG_INFO("phi_hatta", "FF alphaP parameter" + std::to_string(alphaP));
    ff_func_ = [=, *this](double Q2, double W, double t, double Mt) -> double {
      return physics::exp_ff_normalized(Q2, W, t, Mt, vm_.mass(), B0, alphaP);
    };
  } else if (ff == "dipole") {
    const double Mg2 = cf.get<double>(path / "ff" / "Mg2");
    LOG_INFO("phi_hatta", "FF Mg^2 parameter" + std::to_string(Mg2));
    ff_func_ = [=, *this](double Q2, double W, double t, double Mt) -> double {
      return physics::dipole_ff_normalized(Q2, W, t, Mt, vm_.mass(), Mg2);
    };
  } else {
    throw cf.value_error("ff/function", ff);
  }
  LOG_INFO("phi_hatta", "VM: " + std::string(vm_.pdg()->GetName()));
  LOG_INFO("phi_hatta", "recoil: " + std::string(recoil_.pdg()->GetName()));
}

lA_event phi_hatta::generate(const lA_data& initial) {

  // generate a mass() in case of non-zero width, initialize the particles
  particle vm = {vm_.type(), rng()};
  particle recoil = {recoil_.type(), rng()};
  //LOG_INFO("phi_hatta", "HK WAS HERE ");
  // shortcuts
  const auto& gamma = initial.photon();
  const auto& target = initial.target();

  // check if enough energy available
  if (gamma.W2() < threshold2(vm, recoil)) {
    LOG_JUNK("phi_hatta", "Not enough phase space available - W2: " +
                               std::to_string(gamma.W2()) + " < " +
                               std::to_string(threshold2(vm, recoil)));
    return lA_event{0.};
  }

  // generate a phase space point
  const double t = rng()->Uniform(max_t_range_.min, max_t_range_.max);

  LOG_JUNK("phi_hatta", "t: " + std::to_string(t));

  // check if kinematically allowed
  if (physics::t_range(gamma.W2(), gamma.Q2(), target.particle().mass(),
                       vm.mass(), recoil.mass())
          .excludes(t)) {
    LOG_JUNK("phi_hatta", "t outside of the allowed range for this W2")
    return lA_event{0.};
  }

  // evaluate the cross section
  const double R = physics::R_phi_hatta(gamma.Q2(), vm_.mass(), c_R_);
  const double sigmaT =
      physics::sigmaT_phi_hatta(gamma.Q2(), gamma.W(), target.particle().mass(),
                               vm_.mass(), alpha_1_, alpha_2_, alpha_3_, nu_T_);
  const double ff =
      ff_func_(gamma.Q2(), gamma.W(), t, target.particle().mass());
  const double xs = (1 + gamma.epsilon() * R) * sigmaT * ff;

  LOG_JUNK("phi_hatta",
           "xsec: " + std::to_string(xs) + " < " + std::to_string(max_));
  LOG_JUNK("phi_hatta", "sigmaT: " + std::to_string(sigmaT));
  LOG_JUNK("phi_hatta", "R: " + std::to_string(R));
  LOG_JUNK("phi_hatta", "ff: " + std::to_string(ff));

  // return a new VM event
  return make_event(initial, t, vm, recoil, xs, R);
}

// =============================================================================
// phi_hatta::calc_max_xsec(cf)
//
// Utility function for the generator initialization
//
// max cross section as defined by the beam/target and phase-space settings
//
// The cross section is maximum when
//  * W2 is maximum.
//  * no t-requirement (ff smaller than 1 everywhere)
//  * Q2 is zero
// =============================================================================
double phi_hatta::calc_max_xsec(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the
  // lepton beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                        cf.get<double>("beam/lepton/energy")};
  const particle target{initial::estimated_target(cf)};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W
  const double Wmax = opt_W_range ? fmin(opt_W_range->max * opt_W_range->max,
                                         (photon.p() + target.p()).M())
                                  : (photon.p() + target.p()).M();
  return physics::sigmaT_phi_hatta(0.01, Wmax, target.mass() * 1.0001, vm_.mass(),
                                  alpha_1_, alpha_2_, alpha_3_, nu_T_);
} // namespace lA

// =============================================================================
// phi_hatta::calc_max_t_range(cf)
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
interval<double> phi_hatta::calc_max_t_range(const configuration& cf,
                                              const string_path& path) const {
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
  // check if we have a user-defined t-range set
  const auto opt_t_range = cf.get_optional_range<double>(path / "t_range");
  const interval<double> tlim =
      opt_t_range ? interval<double>{fmax(tlim1.min, opt_t_range->min),
                                     fmin(tlim2.max, opt_t_range->max)}
                  : interval<double>{tlim1.min, tlim2.max};
  return tlim;
}
// =============================================================================
// phi_hatta::threshold2()
//
// utility function returns the production threshold squared for the chosen
// particles. Important to re-calculate in case of a particle with non-zero
// widht.
// =============================================================================
double phi_hatta::threshold2(const particle& vm,
                              const particle& recoil) const {

  return recoil.mass2() + vm.mass2() + 2 * vm.mass() * recoil.mass();
}

// =============================================================================
// create the lA_event dataf, calculates the final state four-vectors in
// the lab-frame
// =============================================================================
lA_event phi_hatta::make_event(const lA_data& initial, const double t,
                                particle vm, particle X, const double xs,
                                const double R) {
  const auto& gamma = initial.photon();
  const auto& target = initial.target();

  lA_event e{initial, xs, 1., R};

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
  // 2. rotated TRF -> TRF
  vm.rotate_uz(phot.p());
  X.rotate_uz(phot.p());
  // 3. TRF -> lab
  vm.boost(boost_to_trf.Inverse());
  X.boost(boost_to_trf.Inverse());

  // update VM particle status
  vm.update_status(particle::status_code::UNSTABLE_SCHC);

  // set vertex info
  vm.vertex() = phot.vertex();
  X.vertex() = phot.vertex();

  // add to the event
  e.add_leading(vm, e.photon_index(), e.target_index());
  e.add_recoil(X, e.photon_index(), e.target_index());

  // all done!
  return e;
}
  //FACTORY_REGISTER2(lA::generator, lager::lA::phi_hatta, "phi_hatta");
} // namespace lA
} // namespace lager



