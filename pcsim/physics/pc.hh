#ifndef PCSIM_PHYSICS_PC_DEFINED
#define PCSIM_PHYSICS_PC_DEFINED

#include <pcsim/core/configuration.hh>

#include <TMath.h>

#include <cmath>


// cross section definitions
namespace pcsim {
namespace physics {

// Cross section 
// Values are in units of nb/GeV^2
class pc_xsec : public generator {
public:
  pc_xsec(const ptree& settings, const string_path& path)
      : configurable(settings, path), coupling_{conf().get<double>("coupling")} {}

#if 0
  // This function returns the cross-section for a charmed pentaquark based on a
  // spline fit to the cross section (to be given through a configuration file)
  // as a function of W only, based on Qian Wang's paper: Phys.Rev. D92 (2015)
  // 034022.  It ignores any angular considerations for the cross-section.
  // Returned result is in nb.
  double operator()(const double W) const {
    return coupling_ * spline::eval(W);
  }
#endif
  // crude gaussian fit to the interpolated numbers, needs MAJOR
  // refinements!
  double operator() (const double W) const {
    return coupling_ * 10000. * TMath::Gaus(W, 4.450, .02);
  }

private:
  const double coupling_;
};

} // physics
} // pcsim

#endif
