#ifndef PCSIM_PHYSICS_BREMSSTRAHLUNG_DEFINED
#define PCSIM_PHYSICS_BREMSSTRAHLUNG_DEFINED

#include <cmath>
#include <TSpline.h>

// photon beam intensities
namespace pcsim {
namespace physics {

// Bremsstrahlung spectrum, using the approximate formula by Tsai and Whitis,
// SLAC-PUB-184 1966 (eq. 25), valid for thicker targets.
// parameters:
//  * t: radiation length
//  * E0: primary electron beam energy
//  * k: photon energy
// note: the approximation is valid within 15% as long as 0.2 < k/E0 < 1
//       and t < 2 radiation lengths.
inline double bremsstrahlung_intensity_approx(const double t, const double E0,
                                              const double k) {
  double num = std::pow(1.0 - k / E0, 4.0 * t / 3.0) - std::exp(-7.0 * t / 9.0);
  double den = 7.0 / 9.0 + (4.0 / 3.0) * std::log(1.0 - k / E0);
  return (num / (k * den));
}

// Bremsstrahlung spectrum for a 10% r.l. radiator interpolated between the
// resultes of the exact calculation by Tsai and Whitis (SLAC-PUB-184 1966
// (Table I))
inline double bremsstrahlung_intensity_10_param(const double E0,
                                                const double k) {
  static double xv[] = {0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.82,
                        0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1};
  static double yv[] = {1.15277, 1.04815, 0.96490, 0.90130, 0.85591, 0.82716,
                        0.81289, 0.80929, 0.80921, 0.80908, 0.80873, 0.80789,
                        0.80616, 0.80289, 0.79691, 0.78565, 0.76136, 0.64338};
  static const TSpline3 brems10{"brems10", xv, yv, 18};
  return brems10.Eval(k / E0) * 0.1 / k;
}

} // physics
} // pcsim

#endif
