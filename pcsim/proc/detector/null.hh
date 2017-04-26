#ifndef PCSIM_PROC_DETECTOR_NULL_LOADED
#define PCSIM_PROC_DETECTOR_NULL_LOADED

#include <pcsim/core/generator.hh>
#include <pcsim/gen/lp_gamma_event.hh>

namespace pcsim {
namespace detector {

class null : public event_processor<lp_gamma_event> {
public:
  using base_type = event_processor<lp_gamma_event>
  none(const configuration&, const string_path&, std::shared_ptr<TRandom> r)
      : base_type{r} {}

  virtual void process(lp_gamma_event& e);
};

} // namespace detector
} // namespace pcsim

#endif
