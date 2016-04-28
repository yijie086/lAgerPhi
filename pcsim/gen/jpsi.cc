#include "jpsi.hh"
#include <pcsim/physics/decay.hh>
#include <pcsim/physics/pdg.hh>

namespace pcsim {
namespace gen {

jpsi::jpsi(const ptree& settings, const string_path& path,
           std::shared_ptr<TRandom> r)
    : base_type{settings, path, "t-channel J/#Psi Generator", std::move(r)}
    , xsec_{settings, path / "xsec"}
    , me_{physics::PDG_ELECTRON.Mass()}
    , Mjp_{physics::PDG_JPSI.Mass()}
    , Mp_{physics::PDG_PROTON.Mass()}
    , Wjp_{physics::PDG_JPSI_WIDTH}
    , ctheta_min_{-1.}
    , ctheta_max_{1.} {}

jpsi_event jpsi::gen_impl(const photon_beam& photon) {
  jpsi_event ev;

  // get B-W J/Psi mass
  const double Mj = rng()->BreitWigner(Mjp_, Wjp_);

  // do we have enough energy to continue?
  const double Emin = (Mj * Mj + 2 * Mj * Mp_) / (2 * Mp_);
  if (photon.energy < Emin) {
    // not enough energy available!
    ev.good = false;
    return ev;
  }

  // setup the initial state
  ev.beam.SetPxPyPzE(0, 0, photon.energy, photon.energy);
  ev.target.SetPxPyPzE(0, 0, 0, Mp_);
  ev.flux = photon.flux;

  // scattering in the photon-proton CM frame
  const auto cm = (ev.beam + ev.target);
  const auto beta_cm = cm.BoostVector();
  auto beam_cm = ev.beam;
  beam_cm.Boost(-beta_cm);
  auto target_cm = ev.target;
  target_cm.Boost(-beta_cm);
  const double Ep_cm = target_cm.E();
  const double Pp_cm = target_cm.Vect().Mag();
  // set s and W (note, s == W^2 for t-channel processes
  ev.s = cm.M2();
  ev.W = std::sqrt(ev.s);

  // J/Psi and p' CM energy and momentum
  const double Etot_cm = beam_cm.E() + Ep_cm;
  const double Ej_cm =
      (Etot_cm * Etot_cm + Mj * Mj - Mp_ * Mp_) / (2 * Etot_cm);
  const double Epp_cm =
      (Etot_cm * Etot_cm - Mj * Mj + Mp_ * Mp_) / (2 * Etot_cm);
  const double Pj_cm = std::sqrt(Ej_cm * Ej_cm - Mj * Mj);
  const double Ppp_cm = std::sqrt(Epp_cm * Epp_cm - Mp_ * Mp_);

  // get t, tmin and tmax
  ev.tmin = 2 * Mp_ * Mp_ - 2 * (Ep_cm * Epp_cm - Pp_cm * Ppp_cm * ctheta_min_);
  ev.tmax = 2 * Mp_ * Mp_ - 2 * (Ep_cm * Epp_cm - Pp_cm * Ppp_cm * ctheta_max_);
  ev.t = rng()->Uniform(ev.tmin, ev.tmax);

  // get the J/Psi and p' CM 4-vectors
  const double ctheta_cm =
      (ev.t + 2 * Ep_cm * Epp_cm - 2 * Mp_ * Mp_) / (2 * Pp_cm * Ppp_cm);
  const double theta_cm = std::acos(ctheta_cm);
  const double phi_cm = rng()->Uniform(0.0, TMath::TwoPi());
  TVector3 ptemp;
  ptemp.SetMagThetaPhi(Ppp_cm, theta_cm, phi_cm);
  ev.recoil = {ptemp, Epp_cm};
  ptemp.SetMagThetaPhi(Pj_cm, theta_cm, phi_cm);
  ev.jpsi = {-ptemp, Ej_cm};

  // return to the lab frame
  ev.recoil.Boost(beta_cm);
  ev.jpsi.Boost(beta_cm);

  // get the decay leptons
  physics::decay_jpsi_lepton(rng(), ev.jpsi, me_, ev.positron, ev.electron);

  // get the cross section
  ev.xsec = xsec_(ev.s, ev.t, Mj) * (ev.tmax - ev.tmin);
  ev.weight = ev.xsec;

  // all done!
  ev.good = true;
  return ev;
}

} // gen
} // pcsim
