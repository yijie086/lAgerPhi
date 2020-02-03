#include "detector.hh"
#include <liege/proc/detector/null.hh>

namespace liege {
namespace detector {

// initialize the factory
factory<detector, const configuration&, const string_path&,
        std::shared_ptr<TRandom>>
    detector::factory;

// register our generators
//FACTORY_REGISTER(detector, null, "4pi");

} // namespace detector
} // namespace liege
