#ifndef PCSIM_GEN_DECAY_GAMMA_P_EVENT_LOADED
#define PCSIM_GEN_DECAY_GAMMA_P_EVENT_LOADED

#include <pcsim/core/event.hh>
#include <pcsim/core/decay.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/physics/decay.hh>

namespace pcsim {
namespace decay {
// =============================================================================
// DECAY ALL UNSTABLE PARTICLES IN A GAMMA_P EVENT
//
// also handles chained decays
//
// Supported channels:
//  * Pc according to Wang
//  * e+e- decay of VMs
//
// Note: will add event weight to account for branching ratios when not
// simulating the full decay width
// =============================================================================
// TODO update with proper gamma_p_event instead of the main event
class gamma_p_event : public decay_handler<event> {
public:
  gamma_p_event(std::shared_ptr<TRandom> rng) : decay_handler{std::move(rng)} {}
  virtual void decay(event& e);
};

} // namespace decay
} // namespace pcsim

#endif
