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

#ifndef LAGER_GEN_INITIAL_GENERATOR_LOADED
#define LAGER_GEN_INITIAL_GENERATOR_LOADED

#include <lager/core/factory.hh>
#include <lager/core/generator.hh>
#include <lager/gen/initial/data.hh>

// =============================================================================
// main include file for beam generators
// =============================================================================

namespace lager {
namespace initial {

// =============================================================================
// parent template for beam generators
// =============================================================================
template <class Data, class... Input>
class generator : public lager::generator<Data, Input...> {
public:
  using data_type = Data;
  using base_type = lager::generator<Data, Input...>;

  static factory<generator, const configuration&, const string_path&,
                 std::shared_ptr<TRandom>>
      factory_instance;

  generator(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}
};

template <class Data, class... Input>
factory<generator<Data, Input...>, const configuration&, const string_path&,
        std::shared_ptr<TRandom>>
    generator<Data, Input...>::factory_instance;

// =============================================================================
// beam generator types
// =============================================================================
using vertex_generator = generator<vertex>;
using beam_generator = generator<beam, vertex>;
using target_generator = generator<target, beam>;
using photon_generator = generator<photon, beam, target>;

} // namespace initial
} // namespace lager

#endif
