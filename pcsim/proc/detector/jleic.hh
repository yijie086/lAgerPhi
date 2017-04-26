#ifndef PCSIM_PROC_DETECTOR_JLEIC_LOADED
#define PCSIM_PROC_DETECTOR_JLEIC_LOADED

#include <pcsim/core/generator.hh>
#include <pcsim/gen/lp_gamma_event.hh>
#include <pcsim/proc/detector/eic-fastmc/acceptance.h>
#include <pcsim/proc/detector/eic-fastmc/resolution.h>

namespace pcsim {
namespace detector {

class jleic : public event_processor<lp_gamma_event> {
public:
  using base_type = event_processor<lp_gamma_event>

  jleic(const configuration&, const string_path&, std::shared_ptr<TRandom> r)
      : base_type{r} {}

  virtual void process(lp_gamma_event& e);

private:
  jleic_impl::acceptance acceptance_;
  jleic_impl::resolution resolution_;
};

} // namespace detector
} // namespace pcsim

#endif
