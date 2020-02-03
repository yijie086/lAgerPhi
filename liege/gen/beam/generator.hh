#ifndef LIEGE_GEN_BEAM_GENERATOR_LOADED
#define LIEGE_GEN_BEAM_GENERATOR_LOADED

#include <liege/core/factory.hh>
#include <liege/core/generator.hh>
#include <liege/gen/beam/photon.hh>
#include <liege/gen/beam/primary.hh>

// =============================================================================
// main include file for beam generators
// =============================================================================

namespace liege {
namespace beam {

// =============================================================================
// parent template for beam generators
// =============================================================================
template <class Data, class... Input>
class generator : public liege::generator<Data, Input...> {
public:
  using data_type = Data;
  using base_type = liege::generator<Data, Input...>;

  static factory<generator, const configuration&, const string_path&,
                 std::shared_ptr<TRandom>>
      factory;

  generator(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}
};

template <class Data, class... Input>
factory<generator<Data, Input...>, const configuration&, const string_path&,
        std::shared_ptr<TRandom>>
    generator<Data, Input...>::factory;

// =============================================================================
// beam generator types
// =============================================================================
using primary_generator = generator<primary, particle::XYZTVector>;
using photon_generator = generator<photon, primary, primary>;
using vertex_generator = generator<particle::XYZTVector>;

} // namespace beam
} // namespace liege

#endif
