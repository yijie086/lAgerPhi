#ifndef PCSIM_PROC_DETECTOR_DETECTOR_LOADED
#define PCSIM_PROC_DETECTOR_DETECTOR_LOADED

#include <pcsim/core/event.hh>
#include <pcsim/core/generator.hh>

// =============================================================================
// main include file for detectors
// =============================================================================

namespace pcsim {
namespace detector {

// =============================================================================
// parent class for detectors
// =============================================================================
class detector : public event_processor<event> {
public: 
  using event_type = event;
  using base_type = event_processor<event>;

  static factory<detector, const configuration&, const string_path&,
                 std::shared_ptr<TRandom>>
      factory;

  detector(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}
};

} // namespace detector
} // namespace pcsim

#endif
