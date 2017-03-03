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
    , y_acc_{conf().get<double>("acceptance/y") / 1000.}
    , p_smear_{conf().get<double>("smearing/momentum")}
    , x_smear_{conf().get<double>("smearing/x") / 1000.}
    , y_smear_{conf().get<double>("smearing/y") / 1000.} {
  LOG_INFO(path.str(),
           "Momentum [GeV]: " + conf().get<std::string>("momentum/center"));
  LOG_INFO(path.str(), "Theta [deg.]: " + conf().get<std::string>("theta"));
  LOG_INFO(path.str(), "Charge [e]: " + conf().get<std::string>("charge"));
  LOG_INFO(path.str(), "Momentum Window [GeV]: [" +
                           std::to_string(p_range_.min) + ", " +
                           std::to_string(p_range_.max) + "]");
  LOG_INFO(path.str(), "Acceptance (in, out) [mrad]: [" +
                           conf().get<std::string>("acceptance/x") + ", " +
                           conf().get<std::string>("acceptance/y") + "]");
  if (p_smear_ > 0) {
    LOG_INFO(path.str(), "Momentum Smearing: " +
                             conf().get<std::string>("smearing/momentum"));
  }
  if (x_smear_ > 0 || y_smear_ > 0) {
    LOG_INFO(path.str(), "Angular Smearing (in, out) [mrad]: [" +
                             conf().get<std::string>("smearing/x") + ", " +
                             conf().get<std::string>("smearing/y") + "]");
  }
}

} // gen
} // pcsim
