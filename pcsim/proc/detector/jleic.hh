#ifndef PCSIM_PROC_DETECTOR_JLEIC_LOADED
#define PCSIM_PROC_DETECTOR_JLEIC_LOADED

#include <pcsim/proc/detector/detector.hh>
#include <pcsim/proc/detector/eic-fastmc/acceptance.h>
#include <pcsim/proc/detector/eic-fastmc/resolution.h>

namespace pcsim {
namespace detector {

class jleic : public detector {
public:
  using base_type = detector;

  jleic(const configuration&, const string_path&, std::shared_ptr<TRandom> r)
      : base_type{r}, resolution_{r} {}

  virtual void process(event& e) const;

private:
  mutable jleic_impl::acceptance acceptance_;
  mutable jleic_impl::resolution resolution_;
};

} // namespace detector
} // namespace pcsim

#endif
