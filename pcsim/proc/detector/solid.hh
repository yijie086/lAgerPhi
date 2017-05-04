#ifndef PCSIM_PROC_DETECTOR_SOLID_LOADED
#define PCSIM_PROC_DETECTOR_SOLID_LOADED

#include <pcsim/proc/detector/detector.hh>
#include <memory>
#include <TH2F.h>
#include <TFile.h>

namespace pcsim {
namespace detector {

class solid : public detector {
public:
  using base_type = detector;

  solid(const configuration&, const string_path&, std::shared_ptr<TRandom> r);

  virtual void process(event& e) const;

private:
  std::shared_ptr<TFile> acceptance_file_;
  TH2F* forward_angle_;
  TH2F* large_angle_;
};

} // namespace detector
} // namespace pcsim

#endif
