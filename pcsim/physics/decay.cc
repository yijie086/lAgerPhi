#include "decay.hh"
#include <cmath>

namespace pcsim {
namespace decay {

// two body decay of a particle 'part' into two particles xx (xx.first,
// xx.second), with angles of the first decay particle ('theta1', 'phi1')
//
// Note: the angles are assumed to be in the helicity frame of 'part'
std::pair<particle, particle> decay_2body(const TLorentzVector& part,
                                          const double theta_1,
                                          const double phi_1,
                                          std::pair<particle, particle>& xx) {

  // calculate the CM kinematics
  const double E = part.mass();
  const double M2_1 = xx.first.mass() * xx.first.mass();
  const double M2_2 = xx.second.mass() * xx.second.mass();
  const double E_1 = (E * E + M2_1 - M2_2) / (2 * E);
  const double E_2 = (E * E - M2_1 + M2_2) / (2 * E);
  const double P_1 = sqrt(E_1 * E_1 - M2_1);
  const double P_2 = sqrt(E_2 * E_2 - M2_2);

  // create the 4-vectors in the part helicity frame
  TVector3 mom;
  mom.SetMagThetaPhi(P_1, theta_1, phi_1);
  xx.first.p() = {mom, E_1};
  // decay particle 1 and two are back-to-back in the CM frame
  mom.SetMagThetaPhi(P_2, theta_1, phi_1);
  xx.second.p() = {-mom, E_2};

  // boost back to the rotated version of the original frame pointing in the
  // direction of part
  const particle::Boost b{
      -(particle::XYZTVector{0, 0, part.p().Vect().Mag(), part.p().E()}
            .BoostToCM())};
  xx.first.boost(b);
  xx.first.boost(b);

  // rotate to the original frame
  xx.first.rotate_uz(part.p().Vect());
  xx.second.rotate_uz(part.p().Vect());

  // that's all
}

} // decay
} // pcsim
