#include "generator.hh"
#include <TRandom.h>
#include <memory>
#include <pcsim/core/configuration.hh>
#include <pcsim/gen/lp_gamma/brodsky_2vmX.hh>
#include <pcsim/gen/lp_gamma/gaussian_1X.hh>

namespace pcsim {

// initialize the factory
//factory<lp_gamma::generator, const configuration&, const string_path&,
//        std::shared_ptr<TRandom>>
//    lp_gamma::generator::factory;
namespace lp_gamma {

// register our generators
//FACTORY_REGISTER(generator, brodsky_2vmX, "brodsky_2vmX");
//FACTORY_REGISTER(generator, gaussian_qpq, "gaussian_1qpq");

} // namespace lp_gamma
} // namespace pcsim
