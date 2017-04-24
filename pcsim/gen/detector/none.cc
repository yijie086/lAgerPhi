#include "none.hh"

namespace pcsim {
namespace detect {

FACTORY_REGISTER(detector, none, "dummy");

bool none::accepted(const particle& part) {
  return false;
}
detected_particle none::generate(const particle& part) { return {part, -1}; }

} // namespace detect
} // namespace pcsim
