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

#ifndef LAGER_PROC_DETECTOR_CONE_LOADED
#define LAGER_PROC_DETECTOR_CONE_LOADED

#include <TFile.h>
#include <TH2F.h>
#include <memory>
#include <lager/core/interval.hh>
#include <lager/proc/detector/detector.hh>
#include <vector>

namespace lager {
namespace detector {

class cone : public detector {
public:
  using base_type = detector;

  cone(const configuration&, const string_path&, std::shared_ptr<TRandom> r);

  virtual void process(event& e) const;

private:
  ROOT::Math::PxPyPzMVector detected_track(const particle& part) const;

  const std::string name_;       // spectromter name
  const int id_{0};              // cone ID
  const interval<double> theta_; // cone theta range
  const interval<double> p_;     // cone momentum range
  const std::vector<int> pid_;   // PID info for acceptable particle types
  const double acceptance_{1.};  // flat acceptance
  const double p_smear_{0.};     // optional momentum smearing
  const double theta_smear_{0.}; // optional angle smearing
  const double phi_smear_{0.};   //
};

} // namespace detector
} // namespace lager

#endif
