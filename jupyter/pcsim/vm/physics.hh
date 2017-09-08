#ifndef VM_PHYSICS_LOADED
#define VM_PHYSICS_LOADED

#include <Math/Vector4D.h>
#include <TMath.h>
#include <cmath>

namespace vm {

// All masses are given in units of GeV

// MASSES
constexpr const double M_EL = .0005109989461;
constexpr const double M_P = .938272;
constexpr const double M_JPSI = 3.096916;
constexpr const double M_UPSILON = 9.46030;

// MASSES SQUARED
constexpr const double M2_EL = M_EL * M_EL;
constexpr const double M2_P = M_P * M_P;
constexpr const double M2_JPSI = M_JPSI * M_JPSI;
constexpr const double M2_UPSILON = M_UPSILON * M_UPSILON;

// THRESHOLD IN W
constexpr const double THRESHOLD_JPSI = M_P + M_JPSI;
constexpr const double THRESHOLD_UPSILON = M_P + M_UPSILON;

// utility functions
//
// transform from s/W2 -> Egamma
inline double s2k(const double s) { return (s - M2_P) / (2. * M_P); }
// transform from Egamma -> s/W2
inline double k2s(const double k) { return (M2_P + 2. * M_P * k); }
// calculate generalized nu (q . P/M_P) from a given Q2 and W
inline double calc_nu(const double Q2, const double W) {
  return (W * W - M2_P + Q2) / (2. * M_P);
}
// calculate generalized Ebeam (k . P/M_P) from a given beam and target 4-vector
inline double calc_Ebeam(const ROOT::Math::XYZTVector& beam,
                         const ROOT::Math::XYZTVector& target) {
  return (beam).Dot(target)/target.M();
}
// utility includes
// virtual photon polarization
// Ebeam: generalized (invariant) beam energy (k . P) / M_P
inline double epsilon(const double Q2, const double W, const double Ebeam) {
  const double nu = calc_nu(Q2, W);
  const double nu2 = nu * nu;
  const double m2 = M2_EL;
  const double tworhopp =
      (2. * Ebeam - nu) * (2. * Ebeam - nu) / (nu2 + Q2) + 1 - 4 * m2 / Q2;
  return 1 + (4 * m2 / Q2 - 2) / tworhopp;
}

}

#endif
