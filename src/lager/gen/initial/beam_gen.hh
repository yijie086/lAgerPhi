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

#ifndef LAGER_GEN_INITIAL_BEAM_GEN_LOADED
#define LAGER_GEN_INITIAL_BEAM_GEN_LOADED

#include <lager/core/factory.hh>
#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/initial/data.hh>

namespace lager {
namespace initial {

// =============================================================================
// constant_beam: primary beam with constant energy
//
// constant 4-vector beam
// =============================================================================
class constant_beam : public beam_generator {
public:
  constant_beam(const configuration& cf, const string_path& path,
                std::shared_ptr<TRandom> r)
      : beam_generator{std::move(r)}
      , beam_{cf.get<std::string>(path / "particle_type"),
              cf.get_vector3<particle::XYZVector>(path / "dir"),
              cf.get<double>(path / "energy"), particle::status_code::BEAM} {
    LOG_INFO("initial::constant_beam",
             "type: " + std::string(beam_.pdg()->GetName()));
    LOG_INFO("initial::constant_beam",
             "energy [GeV]: " + std::to_string(beam_.energy()));
  }

  virtual beam generate(const vertex& vx) {
    particle ret = beam_;
    ret.vertex() = vx;
    return {ret};
  }
  virtual double max_cross_section() const { return 1.; }
  virtual double phase_space() const { return 1.; }

protected:
  const particle beam_;
};

} // namespace initial
} // namespace lager

#endif
