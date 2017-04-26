#ifndef PCSIM_GEN_BEAM_BEAM_LOADED
#define PCSIM_GEN_BEAM_BEAM_LOADED

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
class beam_generator : public generator<Data, Input...> {
public:
  using event_type = Event;
  using base_type = generator<Data, Input...>;

  static factory<beam_generator, const configuration&, const string_path&,
                 std::shared_ptr<TRandom>>
      factory;

  beam_generator(std::shared_ptr<TRandom> r) : base_type{std::move(r)} {}
};

// =============================================================================
// beam generator types
// =============================================================================
using primary_generator = beam_gen<primary>;
using photon_generator = beam_gen<photon, primary, primary>;

} // namespace beam
} // namespace pcsim

#endif
