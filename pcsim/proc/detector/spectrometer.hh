#ifndef PCSIM_PROC_DETECTOR_SPECTROMETER_LOADED
#define PCSIM_PROC_DETECTOR_SPECTROMETER_LOADED

#include <TFile.h>
#include <TH2F.h>
#include <memory>
#include <pcsim/core/interval.hh>
#include <pcsim/proc/detector/detector.hh>
#include <vector>

namespace pcsim {
namespace detector {

class spectrometer : public detector {
public:
  using base_type = detector;

  spectrometer(const configuration&, const string_path&,
               std::shared_ptr<TRandom> r);

  virtual void process(event& e) const;

private:
  std::tuple<double, double, double> track_th_in_out_pz(const particle&) const;
  ROOT::Math::PxPyPzMVector detected_track(const particle&, double,
                                           double) const;

  const std::string name_;        // spectromter name
  const int id_{0};               // spectrometer ID
  const double theta0_;           // polar angle of central ray (in rad)
  const double phi0_;             // azimuthal angle of central ray (in rad)
  const double p0_;               // central momentum (in GeV)
  const interval<double> th_in_;  // in-plane angle with central ray (in rad)
  const interval<double> th_out_; // out-of-plane with to central ray (in rad)
  const interval<double> dp_;     // around central momentum for accetable track
  const interval<double> p_;      // momentum range for acceptable track
  const std::vector<int> pid_;    // PID info for acceptable particle types
  const double acceptance_{1.};   // flat acceptance
  const double p_smear_{0.};      // optional momentum smearing
  const double th_in_smear_{0.};  // optional inbending angle smearing
  const double th_out_smear_{0.}; // optional outbending angle smearing
};

} // namespace detector
} // namespace pcsim

#endif
