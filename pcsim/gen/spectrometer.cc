#include "spectrometer.hh"
namespace pcsim {
namespace gen {

spectrometer::spectrometer(const ptree& settings, const string_path& path,
                           std::shared_ptr<TRandom> r)
    : base_type{settings, path, path.str(), r}
    , p_range_{conf().get<double>("momentum/center") *
                   (1 - conf().get<double>("momentum/dmin")),
               conf().get<double>("momentum/center") *
                   (1 + conf().get<double>("momentum/dmax"))}
    , theta_{conf().get<double>("theta") / 180. * TMath::Pi()}
    , x_acc_{conf().get<double>("acceptance/x") / 1000.}
    , y_acc_{conf().get<double>("acceptance/y") / 1000.} {}

} // gen
} // pcsim
