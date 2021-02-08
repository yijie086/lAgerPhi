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

#ifndef LAGER_GEN_BEAM_VERTEX_GEN_LOADED
#define LAGER_GEN_BEAM_VERTEX_GEN_LOADED

#include <lager/core/generator.hh>
#include <lager/core/interval.hh>
#include <lager/core/particle.hh>
#include <lager/gen/initial/data.hh>
#include <lager/gen/initial/generator.hh>

namespace lager {
namespace initial {

// =============================================================================
// vertex at a specific position (default: 0,0,0)
// =============================================================================
class origin_vertex : public vertex_generator {
public:
  origin_vertex(const configuration& cf, const string_path& path,
                std::shared_ptr<TRandom> r)
      : vertex_generator{std::move(r)}
      , vertex_{0, 0, cf.get<double>(path / "vz", 0), 0} {
    LOG_INFO("initial::origin_vertex", "Vertex generator initialized");
  }

  virtual particle::XYZTVector generate() { return vertex_; }
  virtual double max_cross_section() const { return 1.; }
  virtual double phase_space() const { return 1.; }

protected:
  const particle::XYZTVector vertex_;
};

class linear_vertex : public vertex_generator {
public:
  linear_vertex(const configuration& cf, const string_path& path,
                std::shared_ptr<TRandom> r)
      : vertex_generator{std::move(r)}
      , range_{cf.get_range<double>(path / "range")} {
    LOG_INFO("initial::linear_vertex", "Vertex generator initialized");
    LOG_INFO("initial::linear_vertex",
             "Vertex range [cm]: [" + std::to_string(range_.min) + ", " +
                 std::to_string(range_.max) + "]");
  }

  virtual particle::XYZTVector generate() {
    const double vz = rng()->Uniform(range_.min, range_.max);
    LOG_JUNK2("linear_vertex", "Vertex position [cm]: " + std::to_string(vz));
    return {0., 0., vz, 0};
  }
  virtual double max_cross_section() const { return 1.; }
  virtual double phase_space() const { return 1.; }

protected:
  const interval<double> range_;
};

} // namespace initial
} // namespace lager

#endif
