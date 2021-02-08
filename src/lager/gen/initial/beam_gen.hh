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

#ifndef LAGER_GEN_INITIAL_BEAM_GEN_LOADED
#define LAGER_GEN_INITIAL_BEAM_GEN_LOADED

#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/initial/data.hh>
#include <lager/gen/initial/generator.hh>

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
                std::shared_ptr<TRandom> r);
  virtual beam generate(const vertex& vx) {
    particle ret = beam_;
    ret.vertex() = vx;
    return {ret};
  }
  virtual double max_cross_section() const { return 1.; }
  virtual double phase_space() const { return 1.; }

protected:
  particle get_beam(const configuration& cf, const string_path& path) const;
  const particle beam_;
};

} // namespace initial
} // namespace lager

#endif
