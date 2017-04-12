#ifndef PCSIM_PHYSICS_KINEMATICS_LOADED
#define PCSIM_PHYSICS_KINEMATICS_LOADED

#include <cmath>
#include <pcsim/core/interval.hh>
#include <pcsim/core/particle.hh>

namespace pcsim {
namespace physics {

// calculate the allowed Q2 range for a given lepton beam, target and y
//
// The lower bound is given by the kinematics of t-channel scattering,
// while the upper bound is given by either t-channel kinematics, or the
// requirement that the final W be larger than the target mass
interval<double> Q2_range(const particle& beam, const particle& target,
                          const double y) {
  const double E = beam.mom * target.mom / target.mass;
  const double E2 = E * E;
  const double m2 = beam.mass * beam.mass;
  // lower and upper bound from t-channel process on electron leg
  const double comp1 = -2. * m2 + 2. * E2 * (1. - y);
  const double comp2 = 2. * sqrt((E2 - m2) * (E2 * (1. - y) * (1. - y) - m2));
  const double Q2_low = comp1 - comp2;
  const double Q2_high1 = comp1 + comp2;
  // alternative upper bound from requirement that final state has at least the
  // invariant mass of the target mass (W2min = target.mass^2, meaning that 2 M nu = Q2)
  const double Q2_high2 = 2 * target.mass * E * y;
  return {Q2_low, fmin(Q2_high1, Q2_high2)};
}

}
}


#endif
