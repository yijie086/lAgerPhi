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

#ifndef LAGER_PROC_DETECTOR_DETECTOR_LOADED
#define LAGER_PROC_DETECTOR_DETECTOR_LOADED

#include <lager/core/event.hh>
#include <lager/core/generator.hh>

// =============================================================================
// main include file for detectors
// =============================================================================

namespace lager {
namespace detector {

// =============================================================================
// parent class for detectors
// =============================================================================
class detector : public event_processor<event> {
public: 
  using event_type = event;
  using base_type = event_processor<event>;

  static factory<detector, const configuration&, const string_path&,
                 std::shared_ptr<TRandom>>
      factory_instance;

  detector(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}
};

} // namespace detector
} // namespace lager

#endif
