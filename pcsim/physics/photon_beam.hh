#ifndef PCSIM_PHYSICS_PHOTON_BEAM_DEFINED
#define PCSIM_PHYSICS_PHOTON_BEAM_DEFINED

#include <cmath>

// photon beam intensities
namespace pcsim {

// Bremmsstrahlung spectrum, using the approximate formula by Tsai and Whitis,
// SLAC-PUB-184 1966 (eq. 25)
// parameters:
//  * t: radiation length
//  * E0: primary electron beam energy
//  * k: photon energy
// note: the approximation is valid within 15% as long as 0.2 < k/E0 < 1
double photon_beam_approx(const double t, const double E0, const double k) {
  double num = pow(1.0 - k / E0, 4.0 * t / 3.0) - exp(-7.0 * t / 9.0);
  double den = 7.0 / 9.0 + (4.0 / 3.0) * log(1.0 - k / E0);
  return (num / (k * den));
}

}

#endif
