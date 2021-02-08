// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
// 
// This file is part of lAger.
// 
// lAger is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Shoftware Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// lAger is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with lAger.  If not, see <https://www.gnu.org/licenses/>.
// 

#ifndef LAGER_PHYSICS_PHOTON_LOADED
#define LAGER_PHYSICS_PHOTON_LOADED

#include <TLorentzVector.h>
#include <TMath.h>
#include <cmath>
#include <lager/core/particle.hh>
#include <lager/physics/constants.hh> // for ALPHA

// =============================================================================
// PHOTON and VIRTUAL PHOTON physics
// =============================================================================

namespace lager {
namespace physics {

// =============================================================================
// BREMSSTRAHLUNG
// =============================================================================
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

// =============================================================================
// VIRTUAL PHOTON FLUX
// =============================================================================
//
// =============================================================================
// GAMMA_T
//
// transverse virtual photon flux d(logQ2)d(logy)
// this is the default implementation, rather than d/dQ2dy to avoid running
// into machine precision issues when values of Q2 gets small
//
// Based on 
// Budnev, et. al., PhysRept, 1975, Vol15 (4), pp 181-282 
//    equation (6.8) - (6.10)
//    doi:10.1016/0370-1573(75)90009-5.
inline double gamma_t_log(const double Q2, const double y, const particle& beam,
                          const particle& target) {
  // target-rest-frame lepton energy
  const double E = (beam.p()).Dot(target.p()) / target.mass();
  const double E2 = E * E;
  // beam particle masses
  const double m = beam.mass();
  const double m2 = m * m;
  // other kinematic variables
  const double nu = y * E;
  const double nu2 = nu * nu;
  // density matrix element
  const double tworhopp =
      (2. * E - nu) * (2. * E - nu) / (nu2 + Q2) + 1 - 4 * m2 / Q2;
  // jacobian for dnu -> E dy
  // and
  // jacobian for dy -> y d(logy)
  const double jacobian = E * y;
  // putting all together
  const double gamma_t = ALPHA / 2. / TMath::TwoPi() * sqrt(nu2 + Q2) /
                         (E2 - m2) * tworhopp * jacobian;
  return gamma_t;
}
// same as gamma_t_log, but differential in y and Q2
inline double gamma_t(const double Q2, const double y, const particle& beam,
                      const particle& target) {
  return gamma_t_log(Q2, y, beam, target) / Q2 / y;
}

// =============================================================================
// EPSILON
//
// epsilon = Gamma_L / Gamma_T = 2rho++ / rho00
// Based on 
// Budnev, et. al., PhysRept, 1975, Vol15 (4), pp 181-282 
//    equation (6.8) - (6.10)
//    doi:10.1016/0370-1573(75)90009-5.
inline double epsilon(const double Q2, const double y, const particle& beam,
                      const particle& target) {
  // target-rest-frame lepton energy
  const double E = (beam.p()).Dot(target.p()) / target.mass();
  // beam particle masses
  const double m = beam.mass();
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


} // physics
} // lager

#endif
