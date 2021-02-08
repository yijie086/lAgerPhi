// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
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

// Generic fixed target implementation, to be used for photon generators
// that need to be aware of the %RL traversed.

#include <lager/core/configuration.hh>
#include <lager/core/interval.hh>
#include <lager/core/logger.hh>

namespace lager::initial {
class realistic_target {
public:
  realistic_target(const configuration& cf, const string_path& path);
  // Get the total RL based on the vertex position
  double total_rl(const double vz) const {
    double rl = rl_pre_;
    if (target_range_.includes(vz)) {
      rl += rl_per_cm_ * (vz - target_range_.min);
    } else if (vz >= target_range_.max) {
      rl += rl_per_cm_ * target_range_.width();
    }
    LOG_JUNK2("realistic_target", "Calculated RL for z-vertex position at " +
                                      std::to_string(vz) +
                                      " cm [%]: " + std::to_string(rl * 100.));
    return rl;
  }
  double front() const { return target_range_.min; }
  double back() const { return target_range_.max; }
  double length() const { return target_range_.width(); }

private:
  double calc_rl_pre(const configuration& cf, const string_path& path) const;

  const double rl_pre_;                 // radiator+window preceeding the target
  const double rl_per_cm_;              // rl/cm for the target
  const interval<double> target_range_; // target range in cm
};
} // namespace lager::initial
