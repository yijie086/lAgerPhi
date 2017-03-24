#ifndef PCSIM_PHYSICS_DECAY_LOADED
#define PCSIM_PHYSICS_DECAY_LOADED

#include <pcsim/core/particle.h>
#include <pcsim/core/pdg.h>
#include <tuple>

namespace pcsim {
namespace physics {

// =============================================================================
// GENERIC TWO BODY DECAY
// =============================================================================

// two body decay of a particle 'part' into two particles xx (xx.first,
// xx.second), with angles of the first decay particle ('theta1', 'phi1')
//
// Note: the angles are assumed to be in the helicity frame of 'part'
void decay_2body(const particle& part, const double theta_1, const double phi_1,
                 std::pair<particle, particle>& xx)

} // physics
} // pcsim

#endif
