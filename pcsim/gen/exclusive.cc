#include "exclusive.hh"

#define THROW_EXP_BT

namespace pcsim {
namespace gen {

brodsky_tchannel::brodsky_tchannel(const configuration& cf,
                                   const string_path& path,
                                   std::shared_ptr<TRandom> r)
    : exclusive{cf, path, std::move(r)}
    , recoil_{static_cast<pdg_id>(conf().get<int>("recoil_type"))}
    , vm_{static_cast<pdg_id>(conf().get<int>("vm_type"))}
    , threshold2_{recoil_.mass * recoil_.mass + vm_.mass * vm_.mass +
                  2 * vm_.mass * recoil_.mass}
    , b_{conf().get<double>("b")}
    , c2g_{conf().get<double>("c2g")}
    , max_t_range_{calc_max_t_range(cf)}
    , max_exp_bt_range_{exp(b_ * max_t_range_.min), exp(b_ * max_t_range_.max)}
    , max_{calc_max_xsec(cf)} {}

exclusive_data brodsky_tchannel::generate(const photon_data& phdat, const particle& target) {
  exclusive_data event;

  if (phdat.W2 <= threshold2_) {
    event.cross_section = 0;
    return event;
  }

  const auto exp_bt_lim = exp_bt_range(phdat.W2, phdat.Q2, target.mass);

  // check if our t limits are physical (enough energy available, needed
  // because due to rounding we can sometimes end up just below the lower
  // bound)
  if (isnan(exp_bt_lim.min) || isnan(exp_bt_lim.max) ||
      exp_bt_lim.max <= exp_bt_lim.min) {
    // reject the event if not
    event.cross_section = 0;
    return event;
  }
  // throw flat in exp_bt = exp(bt)
#ifdef THROW_EXP_BT
  event.t = std::log(rng()->Uniform(exp_bt_lim.min, exp_bt_lim.max)) / b_;
#else
  const auto t_lim = t_range(phdat.W2, phdat.Q2, target.mass);
  event.t = rng()->Uniform(t_lim.min, t_lim.max);
#endif
 
  // initiate our event with the cross section and phase space
#ifdef THROW_EXP_BT
  event.cross_section = dsigma_dexp_bt(phdat.W2, target.mass) *
                        exp_bt_lim.width() / max_exp_bt_range_.width();
  event.phase_space = max_exp_bt_range_.width();
#else
  event.cross_section = dsigma_dt(phdat.W2, event.t, target.mass) *
                        t_lim.width() / max_t_range_.width();
  event.phase_space = max_t_range_.width();
#endif

  // create our final state particles in the CM frame
  // start with some convenient mass squares
  const double Mt2 = target.mass * target.mass;
  const double Mr2 = recoil_.mass * recoil_.mass;
  const double Mv2 = vm_.mass * vm_.mass;
  // some shortcuts
  const double W2 = phdat.W2;
  const double W = sqrt(W2);
  const double Q2 = phdat.Q2;
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
      (event.t + 2 * Et_cm * Er_cm - Mt2 - Mr2) / (2 * Pt_cm * Pr_cm);
  const double theta_cm = std::acos(ctheta_cm);
  const double phi_cm = rng()->Uniform(0.0, TMath::TwoPi());
  // create the CM 4-vectors
  TVector3 ptemp;
  ptemp.SetMagThetaPhi(Pv_cm, theta_cm, phi_cm);
  event.vm = {static_cast<pdg_id>(vm_.info), {ptemp, Ev_cm}};
  ptemp.SetMagThetaPhi(Pr_cm, theta_cm + TMath::Pi(), phi_cm);
  event.recoil = {static_cast<pdg_id>(recoil_.info), {ptemp, Er_cm}};

  // calculate some nessery boost and rotation vectors
  TLorentzVector targ = target.mom;
  TLorentzVector phot = phdat.photon.mom;
  // to target rest frame
  const auto beta_targ = target.mom.BoostVector();
  targ.Boost(-beta_targ);
  phot.Boost(-beta_targ);
  // describe photon in a target rest frame where the photon moves along the
  // z-axis
  TVector3 phot_dir = phot.Vect().Unit();
  TLorentzVector phot_z = {0, 0, phot.Vect().Mag(), phot.E()};
  // boost to CM frame from this rotated frame
  const auto beta_cm = (targ + phot_z).BoostVector();
  // ok, that's all we need, now let's transform the VM and recoil to the lab
  // frame
  event.vm.mom.Boost(beta_cm);
  event.recoil.mom.Boost(beta_cm);
  event.vm.mom.RotateUz(phot_dir);
  event.recoil.mom.RotateUz(phot_dir);
  event.vm.mom.Boost(beta_targ);
  event.recoil.mom.Boost(beta_targ);

  // all done!
  //std::cout << "t: " << event.t << " " << photon.W2 << " " << max_ << " "
  //          << event.cross_section << std::endl;

  return event;
}

// diffractive (t-channel) J/Psi production cross section dsigma/dt
// Values are in units of nb/GeV^2, from  Brodsky et. al.,
// Phys.Lett.B498:23-28,2001 (http://arxiv.org/abs/hep-ph/0010343)
//
// original formulas:
//  v = 1.0 / (16.0 * TMath::Pi() * pow((s - Mp * Mp), 2));
//  x = (2.0 * Mp * Mj + Mj * Mj) / (s - Mp * Mp);
//
//  A_2g = 6.499e3 * v * pow(1 - x, 2) * pow((s - Mp * Mp), 2) / (Mj * Mj);
//  A_3g =
//    2.894e3 * v * pow(1 - x, 0) * pow((s - Mp * Mp), 2) / (Mj * Mj * Mj * Mj);
//  ep = exp(t * 1.13);
//  xsec = (A_2g (+ A_3g)) * ep
//
//  with Mp the proton mass, Mj the ccbar mass (J/Psi BW mass)
//
//  constants are from eric's fit
//  NOTE: NOT USED, WE USE the expression in dsigma/d(exp bt) instead!
double brodsky_tchannel::dsigma_dt(const double W2, const double t,
                                   const double Mt) const {
  const double Mt2 = Mt * Mt;
  const double Mv = vm_.mass;
  const double Mv2 = vm_.mass * vm_.mass;
  const double x = (2. * Mt * Mv + Mv2) / (W2 - Mt2);
  //const double ep = exp(t * 1.13);
  const double ep = exp(t * b_);
  static const double v = 1. / (16. * TMath::Pi());
  // 2 gluon term
  //double A2g = 6.499e3 * v * (1 - x) * (1 - x) / Mv2;
  const double A2g = c2g_ * v * (1 - x) * (1 - x) / Mv2;
  // return exp* 2gluon term
  return A2g * ep;
}
double brodsky_tchannel::dsigma_dexp_bt(const double W2,
                                        const double Mt) const {
  const double Mt2 = Mt * Mt;
  const double Mv = vm_.mass;
  const double Mv2 = vm_.mass * vm_.mass;
  const double x = (2. * Mt * Mv + Mv2) / (W2 - Mt2);
  static const double v = 1. / (16. * TMath::Pi());
  // 2 gluon term
  const double A2g = c2g_ * v * (1 - x) * (1 - x) / Mv2;
  // additional jacobian for dt -> d(bt)
  const double jacobian = 1 / b_;
  return A2g * jacobian;
}

// max cross section as defined by the beam/target and phase-space settings
double brodsky_tchannel::calc_max_xsec(const configuration& cf) const {
  const particle beam{static_cast<pdg_id>(cf.get<int>("beam/type")),
                      cf.get_vector3<TVector3>("beam/dir"),
                      cf.get<double>("beam/energy")};
  const particle target{static_cast<pdg_id>(cf.get<int>("target/type")),
                        cf.get_vector3<TVector3>("target/dir"),
                        cf.get<double>("target/energy")};
  // maxium cross section when W2 = s
#ifdef THROW_EXP_BT
  return dsigma_dexp_bt((beam.mom + target.mom).M2(), target.mass) * 1.0001;
#else
  return dsigma_dt((beam.mom + target.mom).M2(), max_t_range_.max,
                   target.mass) *
         1.0001;
#endif
}
// max t-range when we have a photon that carries all of the beam energy
interval<double>
brodsky_tchannel::calc_max_t_range(const configuration& cf) const {
  const particle photon{pdg_id::gamma, cf.get_vector3<TVector3>("beam/dir"),
                        cf.get<double>("beam/energy")};
  const particle target{static_cast<pdg_id>(cf.get<int>("target/type")),
                        cf.get_vector3<TVector3>("target/dir"),
                        cf.get<double>("target/energy")};
  return t_range((photon.mom + target.mom).M2(), 0., target.mass);
}

// calculate the t-range as a function of W2, Q2 and target mass Mt
// note: this is also valid for photo-production, where W2 = s and Q2 = 0
interval<double> brodsky_tchannel::t_range(const double W2, const double Q2,
                                           const double Mt) const {
  const double Mt2 = Mt * Mt;
  const double Mr2 = recoil_.mass * recoil_.mass;
  const double Mv2 = vm_.mass * vm_.mass;
  const double W = sqrt(W2);
  // CM relations for energy and momenta
  const double Et_cm = (W2 + Q2 + Mt2) / (2. * W);
  const double Pt_cm = sqrt(Et_cm * Et_cm - Mt2);
  const double Er_cm = (W2 - Mv2 + Mr2) / (2. * W);
  const double Pr_cm = sqrt(Er_cm * Er_cm - Mr2);
  // t_low when recoil direction changes by 90 degrees compared to target,
  // t_high when the recoil in the target direction
  // note: what we call t_high is also refered to as t_min
  const double t_low = Mt2 + Mr2 - 2. * Et_cm * Er_cm - 2. * Pt_cm * Pr_cm;
  const double t_high = Mt2 + Mr2 - 2. * Et_cm * Er_cm + 2. * Pt_cm * Pr_cm;
  // that's all!
  return {t_low, t_high};
}
interval<double> brodsky_tchannel::exp_bt_range(const double W2,
                                                const double Q2,
                                                const double Mt) const {
  auto tlim = t_range(W2, Q2, Mt);
  return {exp(b_ * tlim.min), exp(b_ * tlim.max)};
}

} // gen
} // pcsim
