#ifndef PCSIM_PROC_DETECTOR_NULL_LOADED
#define PCSIM_PROC_DETECTOR_NULL_LOADED

#include <pcsim/proc/detector/detector.hh>

namespace pcsim {
namespace detector {

class null : public detector {
public:
  using base_type = detector;
  null(const configuration&, const string_path&, std::shared_ptr<TRandom> r)
      : base_type{r} {}

  virtual void process(event& e) const;
};

} // namespace detector
} // namespace pcsim

#endif
