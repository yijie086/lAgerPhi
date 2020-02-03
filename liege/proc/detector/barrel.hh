#ifndef LIEGE_PROC_DETECTOR_BARREL_LOADED
#define LIEGE_PROC_DETECTOR_BARREL_LOADED

#include <TFile.h>
#include <TH2F.h>
#include <memory>
#include <liege/core/interval.hh>
#include <liege/proc/detector/detector.hh>
#include <vector>

namespace liege {
namespace detector {

class barrel : public detector {
public:
  using base_type = detector;

  barrel(const configuration&, const string_path&, std::shared_ptr<TRandom> r);

  virtual void process(event& e) const;

private:
  ROOT::Math::PxPyPzMVector barrel::detected_track(const particle& part) const;

  const std::string name_;       // spectromter name
  const int id_{0};              // barrel ID
  const interval<double> theta_; // barrel theta range
  const interval<double> p_;     // barrel momentum range
  const std::vector<int> pid_;   // PID info for acceptable particle types
  const double acceptance_{1.};  // flat acceptance
  const double p_smear_{0.};     // optional momentum smearing
  const double theta_smear_{0.}; // optional angle smearing
  const double phi_smear_{0.};   //
};

} // namespace detector
} // namespace liege

#endif
