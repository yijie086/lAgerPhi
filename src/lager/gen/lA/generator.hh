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

#ifndef LAGER_GEN_LA_GENERATOR_LOADED
#define LAGER_GEN_LA_GENERATOR_LOADED

#include <lager/core/generator.hh>
#include <lager/gen/lA_event.hh>

// =============================================================================
// main include file for lA generators
// =============================================================================
 
namespace lager {
namespace lA {

// =============================================================================
// lA generator type
// =============================================================================
using generator = lager::process_generator<lA_event, lA_data>;

} // namespace beam
} // namespace lager

#endif
