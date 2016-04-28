#ifndef PCSIM_PHYSICS_JPSI_DEFINED
#define PCSIM_PHYSICS_JPSI_DEFINED

#include <TMath.h>
#include <cmath>
#include <pcsim/core/configuration.hh>
#include <pcsim/physics/pdg.hh>

// cross section definitions
namespace pcsim {
namespace physics {

// diffractive (t-channel) J/Psi production cross section dsigma/dt
// Values are in units of nb/GeV^2, from  Brodsky et. al.,
// Phys.Lett.B498:23-28,2001 (http://arxiv.org/abs/hep-ph/0010343)
//
// original formulas:
//  v = 1.0 / (16.0 * TMath::Pi() * pow((s - Mp * Mp), 2));
//  x = (2.0 * Mp * Mj + Mj * Mj) / (s - Mp * Mp);
//
//  A_2g = 6.499e3 * v * pow(1 - x, 2) * pow((s - Mp * Mp), 2) / (Mj * Mj);
//  A_3g =
//    2.894e3 * v * pow(1 - x, 0) * pow((s - Mp * Mp), 2) / (Mj * Mj * Mj * Mj);
//  ep = exp(t * 1.13);
//  xsec = (A_2g (+ A_3g)) * ep
//
//  with Mp the proton mass, Mj the ccbar mass (J/Psi BW mass)
class jpsi_xsec : public configurable {
public:
  jpsi_xsec(const ptree& settings, const string_path& path)
      : configurable{settings, path}
      , enable_3gluon_{conf().get<bool>("enable_3gluon")}
      , Mp_{PDG_PROTON.Mass()}
      , Mp2_{Mp_ * Mp_}
      , v_{1. / (16. * TMath::Pi())} {}

  double operator()(const double s, const double t, const double Mj) const {
    // kinematic components
    const double Mj2 = Mj * Mj;
    const double x = (2.0 * Mp_ * Mj + Mj2) / (s - Mp2_);
    const double ep = exp(t * 1.13);
    // 2 and (optional) 3 gluon term
    return (A_2g(x, Mj2) + A_3g(x, Mj2)) * ep;
  }

private:
  // 2-gluon term
  double A_2g(const double x, const double Mj2) const {
    return 6.499e3 * v_ * (1 - x) * (1 - x) / Mj2;
  }
  // 3-gluon term
  double A_3g(const double x, const double Mj2) const {
    if (enable_3gluon_) {
      return 2.894e3 * v_ / (Mj2 * Mj2);
    } else {
      return 0.;
    }
  }

  const bool enable_3gluon_;
  const double Mp_;
  const double Mp2_;
  const double v_;
};
}
}

#endif
