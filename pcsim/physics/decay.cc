#include "decay.hh"
#include <TMath.h>
#include <pcsim/core/accept_reject.hh>
#include <pcsim/physics/pdg.hh>

namespace pcsim {
namespace physics {

// J/Psi leptonic decay, assuming helicity conservation (1+cos^2\theta_CM).
// arguments:
//  rng: the random generator to be used
//  jpsi: J/Psi 4-vector
//  ml: lepton mass (electron or muon mass)
//  lplus, lminus: reference to the decay lepton 4-vectors
void decay_jpsi_lepton(std::shared_ptr<TRandom> rng, const TLorentzVector& jpsi,
                       const double ml, TLorentzVector& lplus,
                       TLorentzVector& lminus) {
  // work in the j/psi helicity frame
  const double phi = rng->Uniform(0., TMath::TwoPi());
  const double ctheta = accept_reject_1D(
      rng, {-1, 1}, [](const double ctheta) { return 1 + ctheta * ctheta; },
      2.001);
  const double theta = acos(ctheta);
  const double El = .5 * jpsi.M();
  const double pl = std::sqrt(El * El - ml * ml);
  TVector3 mom;
  mom.SetMagThetaPhi(pl, theta, phi);
  lplus = {mom, El};
  lminus = {-mom, El};
  // boost to a rotated lab-frame pointing in the j/psi direction
  const TLorentzVector jpsi_rot{0, 0, jpsi.Vect().Mag(), jpsi.E()};
  const auto beta = jpsi_rot.BoostVector();
  lplus.Boost(beta);
  lminus.Boost(beta);
  // rotate to the regular lab frame
  TVector3 dir{jpsi.Vect().Unit()};
  lplus.RotateUz(dir);
  lminus.RotateUz(dir);
  // all done!
}
void decay_pc_iso(std::shared_ptr<TRandom> rng, const TLorentzVector& pc,
                  TLorentzVector& proton, TLorentzVector& jpsi) {
  static const double Mp = PDG_PROTON.Mass();
  static const double Mj = PDG_JPSI.Mass();

  // work in the Pc CM frame
  const double phi = rng->Uniform(0., TMath::TwoPi());
  const double ctheta = rng->Uniform(-1, 1);
  const double theta = acos(ctheta);
  const double E_cm = pc.M();
  const double Ep = (E_cm * E_cm + Mp * Mp - Mj * Mj) / (2 * E_cm);
  const double Ej = (E_cm * E_cm - Mp * Mp + Mj * Mj) / (2 * E_cm);
  const double Pp = sqrt(Ep * Ep - Mp * Mp);
  const double Pj = sqrt(Ej * Ej - Mj * Mj);
  TVector3 mom;
  mom.SetMagThetaPhi(Pp, theta, phi);
  proton = {mom, Ep};
  mom.SetMagThetaPhi(Pj, theta, phi);
  jpsi = {-mom, Ej};
  // boost to lab frame
  const auto beta = pc.BoostVector();
  proton.Boost(beta);
  jpsi.Boost(beta);
  // all done!
}

} // physics
} // pcsim
