#ifndef PCSIM_CORE_EVENT_LOADED
#define PCSIM_CORE_EVENT_LOADED

#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <vector>

namespace pcsim {

// base event record carries a list of all particles
// derive from this class for more specialized event records
struct event : generator_data {
  double weight = 1;
  int beam_index;
  int target_index;
  std::vector<mc_particle> part;
};

}

#endif
