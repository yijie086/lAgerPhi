#ifndef LIEGE_PROC_RECONSTRUCTION_RECONSTRUCTION_LOADED
#define LIEGE_PROC_RECONSTRUCTION_RECONSTRUCTION_LOADED

#include <liege/core/generator.hh>

namespace liege {
namespace reconstruction {

template <class Event> using reconstruction = event_processor<Event>;
}
} // namespace liege

#endif
