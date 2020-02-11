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

#include "target_gen.hh"

namespace lager {
namespace initial {

// =======================================================================================
// target beam constructor
// =======================================================================================
primary_target::primary_target(const configuration& cf, const string_path& path,
                               std::shared_ptr<TRandom> r)
    : target_generator{std::move(r)}
    , target_{cf.get<std::string>("beam/ion/particle_type")} {
  LOG_INFO("primary_target", "Using primary ion beam as target");
}
target primary_target::generate(const beam& ion) {
  return target::make_target(ion.particle(), target_);
}

} // namespace initial
} // namespace lager
