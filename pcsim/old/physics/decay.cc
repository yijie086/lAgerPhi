#include "decay.hh"
#include <TMath.h>
#include <pcsim/core/accept_reject.hh>

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

} // physics
} // pcsim
