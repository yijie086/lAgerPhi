// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2020 Sylvester Joosten <sjoosten@anl.gov>
// 
// This file is part of lAger.
// 
// lAger is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Shoftware Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// lAger is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with lAger.  If not, see <https://www.gnu.org/licenses/>.
// 

#ifndef LAGER_PROC_DETECTOR_SPECTROMETER_LOADED
#define LAGER_PROC_DETECTOR_SPECTROMETER_LOADED

#include <TFile.h>
#include <TH2F.h>
#include <memory>
#include <lager/core/interval.hh>
#include <lager/proc/detector/detector.hh>
#include <vector>

namespace lager {
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
} // namespace lager

#endif
