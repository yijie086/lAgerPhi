#ifndef PCSIM_GEN_BEAM_GENERATOR_LOADED
#define PCSIM_GEN_BEAM_GENERATOR_LOADED

#include <pcsim/core/factory.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/gen/beam/photon.hh>
#include <pcsim/gen/beam/primary.hh>

// =============================================================================
// main include file for beam generators
// =============================================================================

namespace pcsim {
namespace beam {

// =============================================================================
// parent template for beam generators
// =============================================================================
template <class Data, class... Input>
class generator : public pcsim::generator<Data, Input...> {
public:
  using data_type = Data;
  using base_type = pcsim::generator<Data, Input...>;

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
} // namespace pcsim

#endif
