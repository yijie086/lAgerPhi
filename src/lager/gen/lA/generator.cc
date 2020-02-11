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

#include "generator.hh"
#include <TRandom.h>
#include <memory>

#include <lager/core/configuration.hh>

namespace lager {

// initialize the factory
//factory<lA::generator, const configuration&, const string_path&,
//        std::shared_ptr<TRandom>>
//    lA::generator::factory;
namespace lA {

// register our generators
//FACTORY_REGISTER(generator, brodsky_2vmX, "brodsky_2vmX");
//FACTORY_REGISTER(generator, gaussian_qpq, "gaussian_1qpq");

} // namespace lA
} // namespace lager
