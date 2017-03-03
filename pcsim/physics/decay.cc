#include "decay.hh"
#include <cmath>

namespace pcsim {
namespace physics {
namespace decay {

// two body decay of a particle 'part' into two particles with masses
// 'mass'.first and 'mass'.second and angles of the first decay particle
// ('theta1', 'phi1')
//
// Note: the angles are assumed to be in the helicity frame of 'part'
std::pair<TLorentzVector, TLorentzVector>
two_body(const TLorentzVector& part, const std::pair<double, double>& mass,
         const double theta_1, const double phi_1) {
  // calculate the CM kinematics
  const double E = part.M();
  const double M2_1 = mass.first * mass.first;
  const double M2_2 = mass.second * mass.second;
  const double E_1 = (E * E + M2_1 - M2_2) / (2 * E);
  const double E_2 = (E * E - M2_1 + M2_2) / (2 * E);
  const double P_1 = sqrt(E_1 * E_1 - M2_1);
  const double P_2 = sqrt(E_2 * E_2 - M2_2);

  // create the 4-vectors in the part helicity frame
  TVector3 mom;
  mom.SetMagThetaPhi(P_1, theta_1, phi_1);
  TLorentzVector p_1{mom, E_1};
  // decay particle 1 and two are back-to-back in the CM frame
  mom.SetMagThetaPhi(P_2, theta_1, phi_1);
  TLorentzVector p_2{-mom, E_2};

  // boost back to the rotated version of the original frame pointing in the
  // direction of part
  const TLorentzVector part_rot{0,0, part.Vect().Mag(), part.E()};
  const TVector3 beta = part_rot.BoostVector();
  p_1.Boost(beta);
  p_2.Boost(beta);

  // rotate to the original frame
  TVector3 dir{part.Vect().Unit()};
  p_1.RotateUz(dir);
  p_2.RotateUz(dir);

  // that's all
  return {p_1, p_2};
}

} // decay
} // physics
} // pcsim
