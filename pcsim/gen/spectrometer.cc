#include "spectrometer.hh"
namespace pcsim {
namespace gen {

spectrometer::spectrometer(const configuration& cf, const string_path& path,
                           std::shared_ptr<TRandom> r)
    : base_type{cf, path, path.str(), r}
    , p_range_{conf().get<double>("momentum/center") *
                   (1 - conf().get<double>("momentum/dmin")),
               conf().get<double>("momentum/center") *
                   (1 + conf().get<double>("momentum/dmax"))}
    , theta_{conf().get<double>("theta") / 180. * TMath::Pi()}
    , charge_{conf().get<int>("charge")}
    , x_acc_{conf().get<double>("acceptance/x") / 1000.}
    , y_acc_{conf().get<double>("acceptance/y") / 1000.} {
  LOG_INFO(path.str(),
           "Momentum [GeV]: " + conf().get<std::string>("momentum/center"));
  LOG_INFO(path.str(), "Theta [deg.]: " + conf().get<std::string>("theta"));
  LOG_INFO(path.str(), "Charge [e]: " + conf().get<std::string>("charge"));
  LOG_INFO(path.str(), "Momentum Window [GeV]: [" +
                           std::to_string(p_range_.min) + ", " +
                           std::to_string(p_range_.max) + "]");
  LOG_INFO(path.str(), "Acceptance (x,y) [mrad]: [" +
                           conf().get<std::string>("acceptance/x") + ", " +
                           conf().get<std::string>("acceptance/y") + "]");
}

} // gen
} // pcsim
