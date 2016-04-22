#ifndef PCSIM_PHYSICS_XSEC_DEFINED
#define PCSIM_PHYSICS_XSEC_DEFINED

#include <pcsim/physics/pdg.hh>
#include <pcsim/util/configuration.hh>
#include <pcsim/util/interval.hh>

// ROOT spline class for interpolation between data points
#include <TSpline3.h>
#include <TMath.h>

#include <cmath>


// cross section definitions
namespace pcsim {

// diffractive (t-channel) J/Psi production cross section dsigma/dt
// Values are in units of nb/GeV^2
class jpsi_xsec : public configurable {
public:
  jpsi_xsec(const ptree& settings, const string_path& path)
      : configurable{settings, path}
      , use_3gluon_{conf().get<bool>("enable_3gluon")}
      , MpMj_{PDG_JPSI.Mass() * PDG_PROTON.Mass()}
      , Mp2_{PDG_PROTON.Mass() * PDG_PROTON.Mass()}
      , Mj2_{PDG_JPSI.Mass() * PDG_JPSI.Mass()}
      , Mj4_{Mj2_ * Mj2_}
      , v_{1.0 / (16.0 * TMath::Pi() * pow((s - Mp2_), 2))}
      , x_{(2.0 * MpMj_ + Mj2_) / (s - Mp2_)} {}

  double operator(const double s, const double t) const {
    const double ep = exp(t * 1.13);
    return (A2g(v_, x_, s) + A3g(v_, x_, s)) * ep;
  }

private:
  // 2-gluon term
  double A_2g(const double v, const double x, const double s) const {
    return 6.499e3 * v * pow(1 - x, 2) * pow((s - Mp2_), 2) / (Mj2_);
  }
  // 3-gluon term
  double A_3g(const double v, const double x, const double s) const {
    if (enable_3gluon_) {
      return 2.894e3 * v * pow(1 - x, 0) * pow((s - Mp2_), 2) / (Mj4_);
    } else {
      return 0.
    }
  }

  const bool enable_3gluon_;
  const double MpMj_;
  const double Mp2_;
  const double Mj2_;
  const double Mj4_;
  const double v_;
  const double x_;
}
}

#endif
