#include "jpsi.hh"
#include <pcsim/core/assert.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/physics/constants.hh>
#include <pcsim/physics/decay.hh>
#include <pcsim/physics/pdg.hh>

namespace pcsim {

// utility functions TODO: move to physics
namespace {
// convert s into photon energy (in the target rest frame)
double s2Egamma(const double s) {
  return (s - physics::M2_PROTON) * physics::ONE_OVER_2M_PROTON;
}
double Egamma2s(const double Egamma) {
  return (2 * physics::M_PROTON * Egamma + physics::M2_PROTON);
}
// threshold energy for photo-production off a proton for a particle with mass M
double E_threshold(const double M) {
  return (M * M + 2 * M * physics::M_PROTON) * physics::ONE_OVER_2M_PROTON;
}
double s_threshold(const double M) { return Egamma2s(E_threshold(M)); }
}

namespace gen {

jpsi::jpsi(const configuration& conf, const string_path& path,
           std::shared_ptr<TRandom> r)
    : base_type{conf, path, "t-channel J/#Psi Generator", std::move(r)}
    , brems_{conf, path / "photon_beam"}
    , xsec_{conf, path / "xsec"}
    , s_range_{s_threshold(physics::M_JPSI), brems_.calc_s_range().max}
    , t_range_{calc_t_range(s2Egamma(s_range_.max))}
    , xsec_max_{calc_max_xsec() * brems_.max()} {
  LOG_INFO("jpsi", "s range [GeV2]: [" + std::to_string(s_range_.min) + ", " +
                       std::to_string(s_range_.max) + "]");
  LOG_INFO("jpsi", "t range [GeV2]: [" + std::to_string(t_range_.min) + ", " +
                       std::to_string(t_range_.max) + "]");
  LOG_INFO("jpsi", "J/Psi t-channel generator initialized");
}

jpsi_event jpsi::gen_impl() {

  jpsi_event ev;

  // generate a new phase space point using the cross section
  double test = 0;
  double xsec = 0; // cross section for the generated phase space point
  double Egamma = 0; // photon energy at this phase space point
  do {
    // record this attempt to generate a new event
    ev.ntrials += 1;

    // generate a phase space point
    ev.s = rng()->Uniform(s_range_.min, s_range_.max);
    ev.W = std::sqrt(ev.s);

    // Setup the initial state
    Egamma = s2Egamma(ev.s);

    // do we have enough energy to continue?
    if (Egamma < E_threshold(physics::M_JPSI)) {
      // reject the event if not
      xsec = -1;
      continue;
    }

    // get the t limits and generate t
    const auto tlim = calc_t_range(Egamma);
    ev.t = rng()->Uniform(tlim.min, tlim.max);


    // get the cross section
    xsec = xsec_(ev.s, ev.t, physics::M_JPSI) * tlim.width() / t_range_.width();

    // apply the corresponding bremsstrahlung cross section
    xsec *= brems_(ev.s);

    // record the phase space volume
    ev.volume = xsec_max_ * t_range_.width() * s_range_.width();

    // accept or reject!
    tassert(xsec <= xsec_max_, "Invalid cross section maximum!");

    // accept or reject!
    test = rng()->Uniform(0, xsec_max_);
  } while (test > xsec);

  // finish setting up the initial state
  ev.beam.SetPxPyPzE(0, 0, Egamma, Egamma);
  ev.target.SetPxPyPzE(0, 0, 0, physics::M_PROTON);

  // scattering in the photon-proton CM frame
  const auto cm = (ev.beam + ev.target);
  const auto beta_cm = cm.BoostVector();
  auto beam_cm = ev.beam;
  beam_cm.Boost(-beta_cm);
  auto target_cm = ev.target;
  target_cm.Boost(-beta_cm);
  const double Ep_cm = target_cm.E();
  const double Pp_cm = target_cm.Vect().Mag();

  // p' CM energy and momentum (needed to calculate t)
  const double Etot_cm = beam_cm.E() + Ep_cm;
  const double Epp_cm = (Etot_cm * Etot_cm - physics::M2_JPSI + physics::M2_PROTON) /
           (2 * Etot_cm);
  const double Ppp_cm = std::sqrt(Epp_cm * Epp_cm - physics::M2_PROTON);
  // J/Psi CM energy and momentum
  const double Ej_cm = (Etot_cm * Etot_cm + physics::M2_JPSI - physics::M2_PROTON) /
          (2 * Etot_cm);
  const double Pj_cm = std::sqrt(Ej_cm * Ej_cm - physics::M2_JPSI);

  // calculate the scattering angle theta
  const double ctheta_cm =
      (ev.t + 2 * Ep_cm * Epp_cm - 2 * physics::M2_PROTON) /
      (2 * Pp_cm * Ppp_cm);

  // get the J/Psi and p' CM 4-vectors
  // NOTE:
  //  * theta is the change in angle from the initial state, don't forget proton
  //    originally was flying backwards (already has theta=pi)
  const double theta_cm = std::acos(ctheta_cm);
  const double phi_cm = rng()->Uniform(0.0, TMath::TwoPi());
  TVector3 ptemp;
  ptemp.SetMagThetaPhi(Pj_cm, theta_cm, phi_cm);
  ev.jpsi = {ptemp, Ej_cm};
  ptemp.SetMagThetaPhi(Ppp_cm, theta_cm + TMath::Pi(), phi_cm);
  ev.recoil = {ptemp, Epp_cm};

  // return to the lab frame
  ev.recoil.Boost(beta_cm);
  ev.jpsi.Boost(beta_cm);

  // get the decay leptons
  physics::decay_jpsi_lepton(rng(), ev.jpsi, physics::M_ELECTRON, ev.positron,
                             ev.electron);

  // all done!
  return ev;
}
  
double jpsi::calc_max_xsec() const {
  return xsec_(s_range_.max, t_range_.max, physics::M_JPSI);
}
interval<double> jpsi::calc_t_range(double Egamma) const {
  // kinematics
  TLorentzVector beam{0, 0, Egamma, Egamma};
  TLorentzVector target{0, 0, 0, physics::M_PROTON};
  // Boost to CM frame
  const auto cm = beam + target;
  const auto beta = cm.BoostVector();
  beam.Boost(-beta);
  target.Boost(-beta);
  // proton scattering kinematics
  const double Etot = beam.E() + target.E();
  const double Ep = target.E();
  const double Pp = target.Vect().Mag();
  const double Epp =
      (Etot * Etot - physics::M2_JPSI + physics::M2_PROTON) / (2 * Etot);
  const double Ppp = std::sqrt(Epp * Epp - physics::M2_PROTON);
  // tmin and tmax
  const double tmin{2 * physics::M2_PROTON - 2 * Ep * Epp - 2 * Pp * Ppp};
  const double tmax{2 * physics::M2_PROTON - 2 * Ep * Epp + 2 * Pp * Ppp};
  return {tmin, tmax};
}

} // gen
} // pcsim
