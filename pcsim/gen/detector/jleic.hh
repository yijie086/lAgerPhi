#ifndef PCSIM_DETECTOR_JLEIC_LOADED
#define PCSIM_DETECTOR_JLEIC_LOADED

#include <pcsim/core/detector.hh>
#include <pcsim/gen/detector/eic-fastmc/acceptance.h>
#include <pcsim/gen/detector/eic-fastmc/resolution.h>

namespace pcsim {
namespace detect {

class jleic : public detector {
public:
  jleic(const configuration&, const string_path&, std::shared_ptr<TRandom> r)
      : detector{r} {}

protected:
  virtual bool accepted(const particle& part);
  virtual detected_particle generate(const particle& part);

private:
  jleic_impl::acceptance acceptance_;
  jleic_impl::resolution resolution_;
};

} // namespace detect
} // namespace pcsim

#endif
