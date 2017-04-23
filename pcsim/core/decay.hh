#ifndef PCSIM_CORE_DECAY_LOADED
#define PCSIM_CORE_DECAY_LOADED

#include <TRandom.h>
#include <memory>
#include <pcsim/core/generator.hh>

namespace pcsim {

// =============================================================================
// Base class for all decays for a certain event type (e.g. gamma_p_event)
//
// A decay handler is a type of generator that decays the decay particles, and
// adds them to the input event.
// =============================================================================
template <class EventType> class decay_handler : public generator<void> {
public:
  decay_handler(std::shared_ptr<TRandom> r) : generator<void>{r} {}

  virtual void decay(EventType& e) = 0;

private:
  virtual void generate() {}
  virtual double max_cross_section() const { return 1.; }
  virtual double phase_space() const { return 1.; }
};

} // namespace pcsim

#endif
