#ifndef PCSIM_GEN_BEAM_PRIMARY_LOADED
#define PCSIM_GEN_BEAM_PRIMARY_LOADED

#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>

namespace pcsim {
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
} // namespace pcsim

#endif
