#include "detector.hh"

namespace pcsim {
namespace detect {
factory<detector, const configuration&, const string_path&,
        std::shared_ptr<TRandom>>
    detector::factory;
} // namespace detect
} // namespace pcsim
