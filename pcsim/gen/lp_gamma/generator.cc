#include "generator.hh"
#include <TRandom.h>
#include <memory>
#include <pcsim/core/configuration.hh>
#include <pcsim/gen/beam/brodsky_2vmX.hh>
#include <pcsim/gen/beam/gaussian_1X.hh>

namespace pcsim {
namespace lp_gamma {

// initialize the factory
factory<generator, const configuration&, const string_path&,
        std::shared_ptr<TRandom>>
    generator::factory;

// register our generators
FACTORY_REGISTER(generator, brodsky_2vmX, "brodsky_2vmX");
FACTORY_REGISTER(generator, gaussian_1X, "gaussian_1X");

} // namespace lp_gamma
} // namespace pcsim
