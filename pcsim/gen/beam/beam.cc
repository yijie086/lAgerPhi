#include "beam.hh"

namespace pcsim {
namespace beam {

factory<primary, const configuration&, const string_path&,
        std::shared_ptr<TRandom>>
    primary::factory;
FACTORY_REGISTER(primary, primary, "primary");

} // gen
} // pcsim
