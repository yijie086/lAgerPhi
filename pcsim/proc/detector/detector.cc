#include "detector.hh"
#include <pcsim/proc/detector/jleic.hh>
#include <pcsim/proc/detector/null.hh>

namespace pcsim {
namespace detector {

// initialize the factory
factory<detector, const configuration&, const string_path&,
        std::shared_ptr<TRandom>>
    detector::factory;

// register our generators
//FACTORY_REGISTER(detector, jleic, "jleic");
//FACTORY_REGISTER(detector, null, "4pi");

} // namespace detector
} // namespace pcsim
