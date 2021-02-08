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

#include "detector.hh"
#include <lager/proc/detector/null.hh>

namespace lager {
namespace detector {

// initialize the factory
factory<detector, const configuration&, const string_path&,
        std::shared_ptr<TRandom>>
    detector::factory_instance;

// register our generators
//FACTORY_REGISTER(detector, null, "4pi");

} // namespace detector
} // namespace lager
