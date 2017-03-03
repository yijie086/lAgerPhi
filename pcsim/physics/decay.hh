#ifndef PCSIM_PHYSICS_DECAY_LOADED
#define PCSIM_PHYSICS_DECAY_LOADED

#include <TLorentzVector.h>
#include <tuple>

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
         const double theta_1, const double phi_1);

} // decay
} // physics
} // pcsim

#endif
