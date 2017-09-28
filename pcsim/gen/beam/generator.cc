#include "generator.hh"
#include <pcsim/gen/beam/photon_gen.hh>
#include <pcsim/gen/beam/primary_gen.hh>

namespace pcsim {
namespace beam {

// initialize the factories
//factory<primary_generator, const configuration&, const string_path&,
//        std::shared_ptr<TRandom>>
//    primary_generator::factory;
//factory<photon_generator, const configuration&, const string_path&,
//        std::shared_ptr<TRandom>>
//    photon_generator::factory;

// register our generators
//FACTORY_REGISTER(primary_generator, beam, "beam");

//FACTORY_REGISTER(photon_generator, bremsstrahlung, "bremsstrahlung");
//FACTORY_REGISTER(photon_generator, vphoton, "vphoton");

} // namespace beam
} // namespace pcsim
