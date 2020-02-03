#ifndef LIEGE_PROC_DECAY_DECAY_LOADED
#define LIEGE_PROC_DECAY_DECAY_LOADED

#include <liege/core/generator.hh>

namespace liege {
namespace decay {

template <class Event> using decay = event_processor<Event>;
}
} // namespace liege

#endif
