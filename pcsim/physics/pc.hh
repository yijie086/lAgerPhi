#ifndef PCSIM_PHYSICS_PC_DEFINED
#define PCSIM_PHYSICS_PC_DEFINED

#include <TF1.h>
#include <TMath.h>
#include <cmath>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/interval.hh>
#include <pcsim/core/logger.hh>

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
  pc_xsec(const configuration& cf, const string_path& path)
      : configurable(cf, path)
      , ampl_{conf().get<double>("ampl")}
      , mean_{conf().get<double>("mean")}
      , sigma_{conf().get<double>("width") / 2.} // sigma is FWHM/2
      , coupling_{conf().get<double>("coupling")}
      , max_{calc_maximum()} {
    LOG_INFO("pc_xsec", "Mean: " + std::to_string(mean_) + " GeV");
    LOG_INFO("pc_xsec", "Width: " + std::to_string(sigma_ * 2) + " GeV");
    LOG_INFO("pc_xsec", "Amplitude: " + std::to_string(ampl_));
    LOG_INFO("pc_xsec", "Coupling: " + std::to_string(coupling_));
  }

  // we throw flat in s, hence the cross section is expressed as a function of
  // s = W^2 with Jacobian dW/ds == 1/(2*W)
  double operator() (const double s) const {
    const double W = sqrt(s);
    const double jacobian = 1 / (2 * W);
    return coupling_ * ampl_ * TMath::Gaus(W, mean_, sigma_) * jacobian;
  }

  double max() const {
    return max_;
  }

private:
  // conservative s-range needed to find cross section maximum
  const interval<double> S_RANGE{1., 100.};

  double calc_maximum() const {
    auto func = [this](double* x, double* p) { return this->operator()(x[0]); };
    TF1 f("pc_xsec", func, S_RANGE.min, S_RANGE.max, 0);
    return f.GetMaximum();
  }

  const double ampl_;     // maximum cross section value
  const double mean_;     // cross section mean
  const double sigma_;    // cross section width
  const double coupling_; // branching ratio of Pc -> J/Psi+p channel
  const double max_;      // cross section maximum
};

} // physics
} // pcsim

#endif
