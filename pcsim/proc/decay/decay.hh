#ifndef PCSIM_PROC_DECAY_DECAY_LOADED
#define PCSIM_PROC_DECAY_DECAY_LOADED

#include <pcsim/core/generator.hh>

namespace pcsim {
namespace decay {

template <class Event> using detector = event_processor<Event>;
}
} // namespace pcsim

#endif
