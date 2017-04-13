#include "gamma_p_2vmX.hh"
#include <pcsim/core/logger.hh>
#include <pcsim/physics/photon.hh>
#include <pcsim/physics/vm.hh>

namespace pcsim {
namespace process {

factory<gamma_p_2vmX> gamma_p_2vmX::factory;

// =============================================================================
// Constructor for gamma_p_2vmX data, calculates the final state four-vectors in
// the lab-frame
// =============================================================================
gamma_p_2vmX_data(const beam::photon_data& photon, const beam::data& target,
                  const double t, const particle& vm1, const particle& X1,
                  const double xs, const double phi, const double R,
                  const double epsilon)
    : generator_data{xs}
    , t_{t}
    , xv_{photon.x() + vm1.mass2() / (2 * target.beam().mass())}
    , Q2plusMv2_{photon.Q2() + vm1.mass2()}
    , R_{R}
    , vm_{vm1}
    , X_{X1} {

  // utility shortcuts
  const double W2 = photon.W2();
  const double W = sqrt(W2);
  const double Q2 = photon.Q2();
  const double y = photon.y();
  const double nu = photon.nu();
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
  const double phi_cm = phi;

  // create the CM 4-vectors
  const particle::Polar3DVector p3_v{Pv_cm, theta_cm, phi_cm};
  vm_.p() = {p3_v.X(), p3_v.Y(), p3_v.Z(), Ev_cm};
  const particle::Polar3DVector p3_r{Pr_cm, theta_cm + TMath::Pi(), phi_cm};
  X_.p() = { p3_r.X(), p3_r.Y(), p3_r.Z(), er_cm;

  // calculate some necessary boost and rotation vectors
  particle targ = target.beam();
  particle phot = photon.beam();
  // lab frame to target rest frame (trf)
  particle::Boost boost_trf{targ.BoostToCM()};
  targ.boost(boost_trf);
  phot.boost(boost_trf);
  // describe photon in a target rest frame where the photon moves along the
  // z-axis
  particle::XYZTVector p_phot_z{0, 0, phot.momentum(), phot.energy()};
  // boost to CM frame from this rotated trf
  particle::boost boost_cm{(targ.p() + p_phot_z).BoostToCM()};

  // go to lab frame
  // 1. CM -> rotated TRF
  vm_.boost(-boost_cm);
  X_.boost(-boost_cm);
  // 2. rotated TRF -> TRF
  vm_.rotate_uz(phot.p());
  X_.rotate_uz(phot.p());
  // 3. TRF -> lab
  vm.boost(-boost_trf);
  x_.boost(-boost_trf);
  // all done!
}

// =============================================================================
// Constructor for process::gamma_p_VmX_brodsky
// =============================================================================
gamma_p_2vmX_brodsky::gamma_p_2vmX_brodsky(const configuration& cf,
                                           const string_path& path,
                                           std::shared_ptr<TRandom> r)
    : gamma_p_2vmX{r}
    , recoil_{static_cast<pdg_id>(cf.get<int>(path / "recoil_type"))}
    , vm_{static_cast<pdg_id>(cf.get<int>(path / "vm_type"))}
    , threshold2_{recoil_.mass * recoil_.mass + vm_.mass * vm_.mass +
                  2 * vm_.mass * recoil_.mass}
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
  LOG_INFO("gamma_p_2vmX_brodsky", "t range [GeV^2]: [" +
                                       std::to_string(max_t_range_.min) + ", " +
                                       std::to_string(max_t_range_.max) + "]");
  LOG_INFO("gamma_p_2vmX_brodsky",
           "b parameter [1/GeV^2]: " + std::to_string(photo_b_));
  LOG_INFO("gamma_p_2vmX_brodsky",
           "2-gluon parameter [1/GeV^2]" + std::to_string(photo_c2g_));
  LOG_INFO("gamma_p_2vmX_brodsky",
           "3-gluon parameter [1/GeV^2]" + std::to_string(photo_c3g_));
  LOG_INFO("gamma_p_2vmX_brodsky",
           "R_vm c-parameter: " + std::to_string(R_vm_c_));
  LOG_INFO("gamma_p_2vmX_brodsky",
           "R_vm n-parameter (power): " + std::to_string(R_vm_n_));
  LOG_INFO("gamma_p_2vmX_brodsky",
           "'Dipole' FF power: " + std::to_string(dipole_n_));
  LOG_INFO("gamma_p_2vmX_brodsky", "VM: " + std::string(vm_.pdg->GetName()));
  LOG_INFO("gamma_p_2vmX_brodsky",
           "recoil: " + std::string(recoil_.pdg->GetName()));
}

gamma_p_2vmX_data gamma_p_2vmX_brodsky::generate(const beam::photon_data photon,
                                                 const beam::data& target) {

  // check if enough energy available
  if (photon.W2() < threshold2_) {
    LOG_JUNK("gamma_p_2vmX_brodsky", "Not enough phase space available - W2: " +
                                         std::to_string(photon.W2()) + " < " +
                                         threshold2_);
    return {0.};
  }

  // generate a phase space point
  const double t = std::log(
      rng()->Uniform(max_exp_bt_range_.min, max_exp_bt_range_.max) / photo_b_);

  LOG_JUNK("gamma_p_2vmX_brodsky", "t: " + std::to_string(t));

  // check if kinematically allowed
  if (physics::t_range(photon.W2, photon.Q2, target.mass(), vm_.mass(),
                       recoil_.mass())
          .excludes(t)) {
    LOG_JUNK("gamma_p_2vmX_brodsky",
             "t outside of the allowed range for this W2")
    return {0.};
  }

  // evaluate the cross section
  const double xs_R = R(photon.Q2());
  const double xs_dipole = dipole(photon.Q2());
  const double xs_photo = dsigma_dexp_bt(Q2, target.mass());
  const double xs = (1 + photon.epsilon() * xs_R) * xs_dipole * xs_photo;

  LOG_JUNK("gamma_p_2vmX_brodsky",
           "xsec: " + std::to_string(xs_photo) + " < " + std::to_string(max_));
  LOG_JUNK("gamma_p_2vmX_brodsky", "R: " + std::to_string(xs_R));
  LOG_JUNK("gamma_p_2vmX_brodksy", "dipole: " + std::to_string(xs_dipole));

  // return a new VM event
  return {
      photon, target, t, vm_, recoil_, xs, rng()->Uniform(0, TMath::TwoPi()),
      xs_R};
}


// =============================================================================
// gamma_p_2vmX_brodksy::calc_max_xsec(cf)
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
gamma_p_2vmX_brodsky::calc_max_xsec(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the lepton
  // beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/dir"),
                        cf.get<double>("beam/energy")};
  const particle target{static_cast<pdg_id>(cf.get<int>("target.type")),
                        cf.get_vector3<particle::XYZVector>("target/dir"),
                        cf.get<double>("target/energy")};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W 
  const double W2max = opt_range ? fmin(opt_range_->max * opt_range_->max,
                                        (photon.mom + target.mom).M2())
                                 : (photon.mom + target.mom).M2();
  return dsigma_dexp_bt(W2max, target.mass()) * 1.0001;
}

// =============================================================================
// gamma_p_2vmX_brodksy::calc_max_t_range(cf)
//
// Utility function for the generator initialization
//
// max t-range occurs for:
//  * photon that carries all of the beam energy (or when we have reached the
//    user-defined maximum value of W max) 
//  * maximum Q2 for the given W (or zero for real photons)
// =============================================================================
gamma_p_2vmX_brodsky::calc_max_t_range(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the lepton
  // beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/dir"),
                        cf.get<double>("beam/energy")};
  const particle target{static_cast<pdg_id>(cf.get<int>("target.type")),
                        cf.get_vector3<particle::XYZVector>("target/dir"),
                        cf.get<double>("target/energy")};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W (where the t-range is the largest)
  const double W2max = opt_range ? fmin(opt_range_->max * opt_range_->max,
                                        (photon.mom + target.mom).M2())
                                 : (photon.mom + target.mom).M2();
  // some shortcuts
  const double M_nu = (photon.p()).Dot(target.p());
  const double Q2max = target.mass2() + 2 * M_nu - W2max;
  // return the corresponding t range
  return t_range(W2max, (Q2max > 1e-10 ? Q2max : 0.), target.mass(), vm_.mass(),
                 recoil_.mass());
}
// =============================================================================
// gamma_p_2vmX_brodksy::exp_bt_range()
//
// Utility function that calculates the exp(bt) range for the t-range
// =============================================================================
interval<double> gamma_p_2vmX_brodksy::exp_bt_range(const double W2,
                                                    const double Q2,
                                                    const double Mt) const {
  auto tlim = t_range(W2, Q2, Mt, vm_.mass(), recoil_.mass());
  return {exp(photo_b_ * tlim.min), exp(photo_b_ * tlim.max)};
}
// =============================================================================
// gamma_p_2vmX_brodksy::dsigma_dexp_bt()
// gamma_p_2vmX_brodksy::R()
// gamma_p_2vmX_brodksy::dipole()
//
// Utility functions to calculate the cross section components
// =============================================================================
gamma_p_2vmX_brodksy::dsigma_dexp_bt(const double W2, const double Mt) const {
  return physics::dsigma_dexp_bt_vm_brodsky(W2, Mt, vm_.mass(), photo_b_,
                                            photo_c2g_, photo_c3g_);
}
gamma_p_2vmX_brodksy::R(const double Q2) const {
  return physics::R_vm_martynov(Q2, vm_.mass(), R_vm_c_, R_vm_n_);
}
gamma_p_2vmX_brodksy::dipole(const double Q2) const {
  return physics::dipole_ff_vm(Q2, vm_.mass(), dipole_n_);
}

} // namespace process
} // namespace pcsim
