#include "gamma_p_2vmX.hh"
#include <TMath.h>
#include <pcsim/core/logger.hh>
#include <pcsim/physics/kinematics.hh>
#include <pcsim/physics/photon.hh>
#include <pcsim/physics/vm.hh>

namespace pcsim {
namespace lp_gamma {

// =============================================================================
// Constructor for lp_gamma::brodsky_2vmX
// =============================================================================
brodsky_2vmX::brodsky_2vmX(const configuration& cf, const string_path& path,
                           std::shared_ptr<TRandom> r)
    : gamma_p_2vmX{r}
    , recoil_{static_cast<pdg_id>(cf.get<int>(path / "recoil_type"))}
    , vm_{static_cast<pdg_id>(cf.get<int>(path / "vm_type"))}
    , photo_b_{cf.get<double>(path / "photo_b")}
    , photo_c2g_{cf.get<double>(path / "photo_c2g")}
    , photo_c3g_{cf.get<double>(path / "photo_c3g")}
    , R_vm_c_{cf.get<double>(path / "R_vm_c")}
    , R_vm_n_{cf.get<double>(path / "R_vm_n")}
    , dipole_n_{cf.get<double>(path / "dipole_n")}
    , max_t_range_{calc_max_t_range(cf)}
    , max_exp_bt_range_{exp(photo_b_ * max_t_range_.min),
                        exp(photo_b_ * max_t_range_.max)}
    , max_{calc_max_xsec(cf)} {
  LOG_INFO("brodsky_2vmX", "t range [GeV^2]: [" +
                               std::to_string(max_t_range_.min) + ", " +
                               std::to_string(max_t_range_.max) + "]");
  LOG_INFO("brodsky_2vmX",
           "b parameter [1/GeV^2]: " + std::to_string(photo_b_));
  LOG_INFO("brodsky_2vmX",
           "2-gluon parameter [1/GeV^2]: " + std::to_string(photo_c2g_));
  LOG_INFO("brodsky_2vmX",
           "3-gluon parameter [1/GeV^2]: " + std::to_string(photo_c3g_));
  LOG_INFO("brodsky_2vmX", "R_vm c-parameter: " + std::to_string(R_vm_c_));
  LOG_INFO("brodsky_2vmX",
           "R_vm n-parameter (power): " + std::to_string(R_vm_n_));
  LOG_INFO("brodsky_2vmX", "'Dipole' FF power: " + std::to_string(dipole_n_));
  LOG_INFO("brodsky_2vmX", "VM: " + std::string(vm_.pdg()->GetName()));
  LOG_INFO("brodsky_2vmX", "recoil: " + std::string(recoil_.pdg()->GetName()));
}

lp_gamma_event brodsky_2vmX::generate(lp_gamma_data& initial) {

  // generate a mass() in case of non-zero width, initialize the particles
  particle vm = {vm_.type(), rng()};
  particle recoil = {recoil_.type(), rng()};

  // shortcuts
  const auto& gamma = initial.photon();
  const auto& target = initial.target();

  // check if enough energy available
  if (gamma.W2() < threshold2(vm, recoil)) {
    LOG_JUNK("brodsky_2vmX", "Not enough phase space available - W2: " +
                                 std::to_string(gamma.W2()) + " < " +
                                 std::to_string(threshold2(vm, recoil)));
    return {0.};
  }

  // generate a phase space point
  const double t =
      std::log(rng()->Uniform(max_exp_bt_range_.min, max_exp_bt_range_.max)) /
      photo_b_;

  LOG_JUNK("brodsky_2vmX", "t: " + std::to_string(t));

  // check if kinematically allowed
  if (physics::t_range(gamma.W2(), gamma.Q2(), target.beam().mass(), vm.mass(),
                       recoil.mass())
          .excludes(t)) {
    LOG_JUNK("brodsky_2vmX", "t outside of the allowed range for this W2")
    return {0.};
  }

  // evaluate the cross section
  const double xs_R = R(gamma.Q2());
  const double xs_dipole = dipole(gamma.Q2());
  const double xs_photo = dsigma_dexp_bt(gamma.W2(), target.beam().mass());
  const double xs = (1 + gamma.epsilon() * xs_R) * xs_dipole * xs_photo;

  LOG_JUNK("brodsky_2vmX",
           "xsec: " + std::to_string(xs_photo) + " < " + std::to_string(max_));
  LOG_JUNK("brodsky_2vmX", "R: " + std::to_string(xs_R));
  LOG_JUNK("brodsky_2vmX", "dipole: " + std::to_string(xs_dipole));

  // return a new VM event
  return make_event(initial, t, vm, recoil, xs, xs_R);
}

// =============================================================================
// brodsky_2vmX::calc_max_xsec(cf)
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
double brodsky_2vmX::calc_max_xsec(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the
  // lepton beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/dir"),
                        cf.get<double>("beam/energy")};
  const particle target{
      static_cast<pdg_id>(cf.get<int>("target/particle_type")),
      cf.get_vector3<particle::XYZVector>("target/dir"),
      cf.get<double>("target/energy")};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W
  const double W2max = opt_W_range ? fmin(opt_W_range->max * opt_W_range->max,
                                          (photon.p() + target.p()).M2())
                                   : (photon.p() + target.p()).M2();
  return dsigma_dexp_bt(W2max, target.mass()) * 1.0001;
}

// =============================================================================
// brodsky_2vmX::calc_max_t_range(cf)
//
// Utility function for the generator initialization
//
// max t-range occurs for:
//  * photon that carries all of the beam energy (or when we have reached the
//    user-defined maximum value of W max)
//  * maximum Q2 for the given W (or zero for real photons)
// =============================================================================
interval<double> brodsky_2vmX::calc_max_t_range(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the
  // lepton beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/dir"),
                        cf.get<double>("beam/energy")};
  const particle target{
      static_cast<pdg_id>(cf.get<int>("target/particle_type")),
      cf.get_vector3<particle::XYZVector>("target/dir"),
      cf.get<double>("target/energy")};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W (where the t-range is the largest)
  const double W2max = opt_W_range ? fmin(opt_W_range->max * opt_W_range->max,
                                          (photon.p() + target.p()).M2())
                                   : (photon.p() + target.p()).M2();
  // some shortcuts
  const double M_nu = (photon.p()).Dot(target.p());
  const double Q2max = target.mass2() + 2 * M_nu - W2max;
  // return the corresponding t range
  // in case of particles with non-zero width, we use M - 4 x sigma
  return physics::t_range(W2max, (Q2max > 1e-10 ? Q2max : 0.), target.mass(),
                          vm_.pole_mass() - vm_.width() * 4.,
                          recoil_.pole_mass() - recoil_.width() * 4);
}
// =============================================================================
// brodsky_2vmX::exp_bt_range()
//
// Utility function that calculates the exp(bt) range for the t-range
// =============================================================================
interval<double> brodsky_2vmX::exp_bt_range(const double W2, const double Q2,
                                            const double Mt) const {
  auto tlim = physics::t_range(W2, Q2, Mt, vm_.mass(), recoil_.mass());
  return {exp(photo_b_ * tlim.min), exp(photo_b_ * tlim.max)};
}
// =============================================================================
// brodsky_2vmX::dsigma_dexp_bt()
// brodsky_2vmX::R()
// brodsky_2vmX::dipole()
//
// Utility functions to calculate the cross section components
// =============================================================================
double brodsky_2vmX::dsigma_dexp_bt(const double W2, const double Mt) const {
  return physics::dsigma_dexp_bt_vm_brodsky(W2, Mt, vm_.mass(), photo_b_,
                                            photo_c2g_, photo_c3g_);
}
double brodsky_2vmX::R(const double Q2) const {
  return physics::R_vm_martynov(Q2, vm_.mass(), R_vm_c_, R_vm_n_);
}
double brodsky_2vmX::dipole(const double Q2) const {
  return physics::dipole_ff_vm(Q2, vm_.mass(), dipole_n_);
}
// =============================================================================
// brodsky_2vmX::threshold2()
//
// utility function returns the production threshold squared for the chosen
// particles. Important to re-calculate in case of a particle with non-zero
// widht.
// =============================================================================
double brodsky_2vmX::threshold2(const particle& vm,
                                const particle& recoil) const {

  return recoil.mass2() + vm.mass2() + 2 * vm.mass() * recoil.mass();
}

// =============================================================================
// create the lp_gamma_event dataf, calculates the final state four-vectors in
// the lab-frame
// =============================================================================
lp_gamma_event brodsky_2vmX::make_event(const lp_gamma_data& initial,
                                        const double t, particle vm1,
                                        particle X1, const double xs,
                                        const double R) {
  const auto& gamma = initial.photon();
  const auto& target = initial.target();

  event e{initial, xs, 1., R};

  // utility shortcuts
  const double W2 = gamma.W2();
  const double W = sqrt(W2);
  const double Q2 = gamma.Q2();
  const double y = gamma.y();
  const double nu = gamma.nu();
  const double Mt2 = target.beam().mass2();
  const double Mr2 = X_.mass2();
  const double Mv2 = vm_.mass2();

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
  particle targ = target.beam();
  particle phot = photon.beam();
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

  // add to the event
  e.add_leading(vm, e.photon_index(), e.target_index());
  e.add_recoil(X, e.photon_index(), e.target_index());

  // all done!
  return e;
}

} // namespace lp_gamma
} // namespace pcsim
