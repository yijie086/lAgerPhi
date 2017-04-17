#ifndef PCSIM_PHYSICS_DECAY_LOADED
#define PCSIM_PHYSICS_DECAY_LOADED

#include <pcsim/core/event.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <tuple>

namespace pcsim {
namespace physics {
// =============================================================================
// DECAY ALL UNSTABLE PARTICLES IN AN EVENT
//
// also handles chained decays
//
// Supported channels:
//  * Pc according to Wang
//  * e+e- decay of VMs
//
// Note: will add event weight to account for branching ratios when not
// simulating the full decay width
//
// TODO make this its own generator
// =============================================================================
class decay_handler : public generator<void>{
public:
  decay_handler(std::shared_ptr<TRandom> rng)
      : generator<void>{std::move(rng)} {}
  // do nothing
  virtual void generate() {}
  virtual double max_cross_section() const {return 1.;}
  virtual double phase_space() const { return 1.; }
  void process(event& e);
};

//(event& e, std::shared_ptr<TRandom> rng);

// =============================================================================
// GENERIC TWO BODY DECAY
// =============================================================================

// two body decay of a particle 'part' into two particles xx (xx.first,
// xx.second), with angles of the first decay particle ('theta1', 'phi1')
//
// Note: the angles are assumed to be in the helicity frame of 'part'
void decay_2body(const particle& part, const double theta_1, const double phi_1,
                 std::pair<particle, particle>& xx);

} // physics
} // pcsim

#endif
