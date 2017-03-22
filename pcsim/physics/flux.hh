#ifndef PCSIM_PHYSICS_FLUX_LOADED
#define PCSIM_PHYSICS_FLUX_LOADED

#include <TLorentzVector.h>
#include <TMath.h>
#include <cmath>
#include <pcsim/physics/constants.hh>

namespace pcsim {
namespace physics {
namespace flux {

// approximate bremsstrahlung d/dk for a radiation length t
// using the approximate formula by Tsai and Whitis,
// SLAC-PUB-184 1966 (eq. 25) valid for thicker targets.
// parameters:
//  * t: radiation length
//  * E0: primary beam energy
//  * k: photon energy
inline double bremsstrahlung_approx(const double t, const double E0, const double k) {
  double num = std::pow(1.0 - k / E0, 4.0 * t / 3.0) - std::exp(-7.0 * t / 9.0);
  double den = 7.0 / 9.0 + (4.0 / 3.0) * std::log(1.0 - k / E0);
  return (num / (k * den));
}

// transverse virtual photon flux d(logQ2)d(logy)
//
// this is the default implementation, rather than d/dQ2dy to avoid running
// into machine precision issues when values of Q2 get very small
inline double gamma_t_log(const double Q2, const double y,
                         const TLorentzVector& beam,
                         const TLorentzVector& target) {
  // target-rest-frame lepton energy
  const double E = beam * target / target.M();
  const double E2 = E * E;
  // beam particle masses
  const double m = beam.M();
  const double m2 = m * m;
  // other kinematic variables
  const double nu = y * E;
  const double nu2 = nu * nu;
  // density matrix element
  const double tworhopp =
      (2. * E - nu) * (2. * E - nu) / (nu2 + Q2) + 1 - 4 * m2 / Q2;
  // jacobian for dnu -> E dy
  // jacobian for dy -> y d(logy)
  const double jacobian = E * y;
  // putting all together
  const double gamma_t = ALPHA / 2. / TMath::TwoPi() * sqrt(nu2 + Q2) /
                         (E2 - m2) * tworhopp * jacobian;
  return gamma_t;
}
// epsilon = Gamma_L / Gamma_T = 2rho++ / rho00
inline double epsilon(const double Q2, const double y,
                      const TLorentzVector& beam,
                      const TLorentzVector& target) {
  // target-rest-frame lepton energy
  const double E = beam * target / target.M();
  // beam particle masses
  const double m = beam.M();
  const double m2 = m * m;
  // other kinematic variables
  const double nu = y * E;
  const double nu2 = nu * nu;
  // density matrix element
  const double tworhopp =
      (2. * E - nu) * (2. * E - nu) / (nu2 + Q2) + 1 - 4 * m2 / Q2;
  // const double rho00 = towrhopp + 4 * m2 / Q2 - 2;
  return 1 + (4 * m2 / Q2 - 2) / tworhopp;
}

// transverse virtual photon flux d/dQ2dy
inline double gamma_t(const double Q2, const double y,
                      const TLorentzVector& beam,
                      const TLorentzVector& target) {
  return gamma_t_log(Q2, y, beam, target) / Q2 / y;
}

} // flux
} // physics
} // pcsim

#endif
