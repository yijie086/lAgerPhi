#ifndef PCSIM_PHYSICS_PC_DEFINED
#define PCSIM_PHYSICS_PC_DEFINED

#include <TMath.h>
#include <cmath>
#include <pcsim/core/configuration.hh>

// cross section definitions
namespace pcsim {
namespace physics {

// Cross section: using a crude gaussian fit to the the results from
//    Wang, Qian, Xiao-Hai Liu, and Qiang Zhao. 2015. “Photoproduction of Hidden
//    Charm Pentaquark States Pc+(4380)and Pc+(4450).” Physical Review D 92 (3):
//    034022–27. doi:10.1103/PhysRevD.92.034022.
// numbers, needs MAJOR refinements!
// Values are in units of nb/GeV^2
class pc_xsec : public configurable {
public:
  pc_xsec(const ptree& settings, const string_path& path)
      : configurable(settings, path)
      , max_{conf().get<double>("max")}
      , mean_{conf().get<double>("mean")}
      , sigma_{conf().get<double>("width") / 2.}
      , coupling_{conf().get<double>("coupling")} {}

  double operator() (const double W) const {
    return coupling_ * max_ * TMath::Gaus(W, mean_, sigma_);
  }

private:
  const double max_;      // maximum cross section value
  const double mean_;     // cross section mean
  const double sigma_;    // cross section width
  const double coupling_; // branching ratio of Pc -> J/Psi+p channel
};

} // physics
} // pcsim

#endif
