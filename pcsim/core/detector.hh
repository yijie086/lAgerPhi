#ifndef PCSIM_CORE_DETECTOR_LOADED
#define PCSIM_CORE_DETECTOR_LOADED

#include <TRandom.h>
#include <memory>
#include <pcsim/core/factory.hh>
#include <pcsim/core/event.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>

namespace pcsim {
namespace detect {

// =============================================================================
// Base class for all detectors
//
// A detector is a type of generator that generates a detected particle from a
// generated particle. The interface is through the detect(event)
// function that will automatically add the detected tracks to the event
// =============================================================================
class detector : public generator<detected_particle, particle> {
public:
  using base_type = generator<detected_particle, particle>;

  static factory<detector, const configuration&, const string_path&,
                 std::shared_ptr<TRandom>>
      factory;

  detector(std::shared_ptr<TRandom> r) : base_type{r} {}

  void detect(event& e) {
    for (const auto& part : e) {
      if (part.final_state() && accepted(part)) {
        e.add_detected(generate(part));
      }
    }
  }

protected:
  virtual detected_particle generate(const particle& part) = 0;
  virtual bool accepted(const particle& part) = 0;
  virtual double max_cross_section() const { return 1.; }
  virtual double phase_space() const { return 1.; }
};

} // namespace detect
} // namespace pcsim

#endif
