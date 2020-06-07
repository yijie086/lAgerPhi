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

#include "beam_gen.hh"

#include <cmath>

namespace lager {
namespace initial {

// =======================================================================================
// constant beam
// =======================================================================================
constant_beam::constant_beam(const configuration& cf, const string_path& path,
                             std::shared_ptr<TRandom> r)
    : beam_generator{std::move(r)}, beam_{get_beam(cf, path)} {
  LOG_INFO("initial::constant_beam",
           "type: " + std::string(beam_.pdg()->GetName()));
  LOG_INFO("initial::constant_beam",
           "energy [GeV]: " + std::to_string(beam_.energy()));
}
particle constant_beam::get_beam(const configuration& cf,
                                 const string_path& path) const {
  const std::string name = cf.get<std::string>(path / "particle_type");
  const auto dir = cf.get_vector3<particle::XYZVector>(path / "dir");
  auto E = cf.get_optional<double>(path / "energy");
  double calc_E = 0;
  if (!E) {
    auto p = cf.get_optional<double>(path / "momentum");
    tassert(E || p, "Either beam energy or momentum should be defined");
    tassert(*p >= 0, "Beam momentum cannot be a negative number");
    particle tmp{name};
    calc_E = std::sqrt(*p * *p + tmp.mass2());
  } else {
    calc_E = *E;
  }
  return {name, dir, calc_E, particle::status_code::BEAM};
}

} // namespace initial
} // namespace lager
