#ifndef LIEGE_PROC_DETECTOR_COMPOSITE_LOADED
#define LIEGE_PROC_DETECTOR_COMPOSITE_LOADED

#include <memory>
#include <liege/core/configuration.hh>
#include <liege/proc/detector/detector.hh>
#include <vector>

namespace liege {
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
} // namespace liege

#endif
