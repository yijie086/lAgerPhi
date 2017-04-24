#ifndef PCSIM_DETECTOR_NONE_LOADED
#define PCSIM_DETECTOR_NONE_LOADED

#include <pcsim/core/detector.hh>

namespace pcsim {
namespace detect {

class none : public detector {
public:
  none(const configuration&, const string_path&, std::shared_ptr<TRandom> r)
      : detector{r} {}

protected:
  virtual bool accepted(const particle& part);
  virtual detected_particle generate(const particle& part);

};

} // namespace detect
} // namespace pcsim

#endif
