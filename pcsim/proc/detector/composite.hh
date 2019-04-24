#ifndef PCSIM_PROC_DETECTOR_COMPOSITE_LOADED
#define PCSIM_PROC_DETECTOR_COMPOSITE_LOADED

#include <memory>
#include <pcsim/core/configuration.hh>
#include <pcsim/proc/detector/detector.hh>
#include <vector>

namespace pcsim {
namespace detector {

class composite : public detector {
public:
  using base_type = detector;

  composite(const configuration&, const string_path&, std::shared_ptr<TRandom> r);

  virtual void process(event& e) const {
    for (const auto& det : detectors_) {
      det->process(e);
    }
  }

private:
  std::vector<std::shared_ptr<detector>> detectors_;
};

} // namespace detector
} // namespace pcsim

#endif
