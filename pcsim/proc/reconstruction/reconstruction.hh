#ifndef PCSIM_PROC_RECONSTRUCTION_RECONSTRUCTION_LOADED
#define PCSIM_PROC_RECONSTRUCTION_RECONSTRUCTION_LOADED

#include <pcsim/core/generator.hh>

namespace pcsim {
namespace reconstruction {

template <class Event> using reconstruction = event_processor<Event>;
}
} // namespace pcsim

#endif
