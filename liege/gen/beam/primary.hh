#ifndef LIEGE_GEN_BEAM_PRIMARY_LOADED
#define LIEGE_GEN_BEAM_PRIMARY_LOADED

#include <liege/core/generator.hh>
#include <liege/core/particle.hh>

namespace liege {
namespace beam {

// =============================================================================
// beam::primary
//
// generic primary beam data
// =============================================================================
class primary : public generator_data {
public:
  primary() = default;
  primary(const particle& part) : beam_{part} {}
  primary(const particle& part, const double xs)
      : generator_data{xs}, beam_{part} {}
  primary(const double xs) : generator_data{xs} {}

  particle& beam() { return beam_; }
  const particle& beam() const { return beam_; }

private:
  particle beam_;
};

} // namespace beam
} // namespace liege

#endif
