#ifndef PCSIM_PHYSICS_BREMSSTRAHLUNG_DEFINED
#define PCSIM_PHYSICS_BREMSSTRAHLUNG_DEFINED

#include <TSpline.h>
#include <cmath>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/interval.hh>
#include <pcsim/physics/constants.hh>

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

// Bremsstrahlung spectrum for a 1%, 5% and 10% r.l. radiator interpolated
// between the  resultes of the exact calculation by Tsai and Whitis
// (SLAC-PUB-184 1966  (Table I))
inline double bremsstrahlung_intensity_001_param(const double E0,
                                                const double k) {
  static double xv[] = {0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.82,
                        0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1};
  static double yv[] = {1.20417, 1.10071, 1.01740, 0.95406, 0.91055, 0.88668,
                        0.88222, 0.89671, 0.90181, 0.90760, 0.91407, 0.92119,
                        0.92889, 0.93710, 0.94566, 0.95421, 0.96164, 0.95451};
  static const TSpline3 brems10{"brems001", xv, yv, 18};
  return brems10.Eval(k / E0) * 0.01 / k;
}
inline double bremsstrahlung_intensity_005_param(const double E0,
                                                const double k) {
  static double xv[] = {0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.82,
                        0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1};
  static double yv[] = {1.18110, 1.07713, 0.99387, 0.93046, 0.88616, 0.86013,
                        0.85123, 0.85735, 0.86001, 0.86302, 0.86625, 0.86956,
                        0.87272, 0.87534, 0.87674, 0.87538, 0.86653, 0.79765};
  static const TSpline3 brems10{"brems005", xv, yv, 18};
  return brems10.Eval(k / E0) * 0.05 / k;
}
inline double bremsstrahlung_intensity_010_param(const double E0,
                                                const double k) {
  static double xv[] = {0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.82,
                        0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1};
  static double yv[] = {1.15277, 1.04815, 0.96490, 0.90130, 0.85591, 0.82716,
                        0.81289, 0.80929, 0.80921, 0.80908, 0.80873, 0.80789,
                        0.80616, 0.80289, 0.79691, 0.78565, 0.76136, 0.64338};
  static const TSpline3 brems10{"brems010", xv, yv, 18};
  return brems10.Eval(k / E0) * 0.1 / k;
}

// Configurable bremsstrahlung spectrum as a function of s, using either the
// approximate formula, the 10% r.l. parameterization, or a flat spectrum
class bremsstrahlung : public configurable {
public:
  // Bremsstrahlung models
  enum class model { FLAT, PARAM_001, PARAM_005, PARAM_010, APPROX };

  bremsstrahlung(const configuration& conf, const string_path& path);

  // get the maximum cross section, useful for accept-reject MC
  double max() const { return max_; }

  // calculate the s-range corresponding to the photon energy range
  interval<double> calc_s_range() const;

  // evaluate the cross section
  double operator()(const double s) const {
    const double Egamma = (s - physics::M_PROTON) * physics::ONE_OVER_2M_PROTON;
    const double jacobian = physics::ONE_OVER_2M_PROTON;
    double xsec = 0;
    if (model_ == model::FLAT) {
      xsec = 1. / jacobian;
    } else if (model_ == model::PARAM_001) {
      xsec = physics::bremsstrahlung_intensity_001_param(E0_, Egamma);
    } else if (model_ == model::PARAM_005) {
      xsec = physics::bremsstrahlung_intensity_005_param(E0_, Egamma);
    } else if (model_ == model::PARAM_010) {
      xsec = physics::bremsstrahlung_intensity_010_param(E0_, Egamma);
    } else {
      xsec = physics::bremsstrahlung_intensity_approx(rl_, E0_, Egamma);
    }
    return xsec * jacobian;
  }

private:
  // utility functions for the constructor
  model get_model() const;
  double calc_max() const;

  // settings
  const model model_; // bremsstrahlung model
  const double rl_;   // number of radiation lenghts (when using approx model,
                      // set to zero otherwise)
  const double E0_;   // primary electron beam energy
  const interval<double> range_; // photon energy range
  // constants
  const double max_;             // the maximum intensity.
};

} // physics
} // pcsim

#endif
