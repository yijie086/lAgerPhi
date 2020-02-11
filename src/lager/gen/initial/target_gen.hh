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

#ifndef LAGER_GEN_INITIAL_TARGET_GEN_LOADED
#define LAGER_GEN_INITIAL_TARGET_GEN_LOADED

#include <TRandom.h>
#include <lager/core/factory.hh>
#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/initial/data.hh>
#include <lager/gen/initial/generator.hh>
#include <memory>

namespace lager {
namespace initial {

// simple target identical to the ion beam
class primary_target : public target_generator {
public:
  primary_target(const configuration& cf, const string_path& path,
                 std::shared_ptr<TRandom> r);

  virtual target generate(const beam& ion);
  // no contribution to the cross section and phase space
  virtual double max_cross_section() const { return 1; }
  virtual double phase_space() const { return 1; }

private:
  const particle target_;
};

} // namespace initial
} // namespace lager

#endif
