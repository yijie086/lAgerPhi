#include "jpsi.hh"
#include <pcsim/core/assert.hh>
#include <pcsim/core/logger.hh>

// cross section definitions
namespace pcsim {
namespace physics {
jpsi_xsec::jpsi_xsec(const configuration& cf, const string_path& path)
    : configurable{cf, path}
    , model_{conf().get<std::string>("type")}
    , enable_3gluon_{model_ == "3gluon"}
    , Mp_{PDG_PROTON.Mass()}
    , Mp2_{Mp_ * Mp_}
    , v_{1. / (16. * TMath::Pi())} {
  LOG_INFO("jpsi_xsec", "Cross section model: " + model_);
  tassert(model_ == "2gluon" || model_ == "3gluon",
          "Invalide J/Psi cross section model");
}

}
}
