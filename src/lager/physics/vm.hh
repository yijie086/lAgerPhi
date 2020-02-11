// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2020 Sylvester Joosten <sjoosten@anl.gov>
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

#ifndef LAGER_PHYSICS_VM_LOADED
#define LAGER_PHYSICS_VM_LOADED

#include <TMath.h>
#include <cmath>

// =============================================================================
// VM Physics routines
// =============================================================================

namespace lager {
namespace physics {

// =============================================================================
// R VM parameterization 
//
// from 
//      Martynov, et. al., “Photoproduction of Vector Mesons in the Soft Dipole
//      Pomeron Model.” PRD 67 (7), 2003. doi:10.1103/PhysRevD.67.074023.
// eq. 31.
//
// J/psi parameters c(or a) and n 
//  * c: 2.164
//  * n: 2.131
// from eq. 18:
//
//      R. Fiore et al., "Exclusive Jpsi electroproduction in a dual model."
//      PRD80:116001, 2009"
//
// Arguments:
//  * Q2: -photon mass squared
//  * Mv: VM mass
//  * c: multiplicative parameter
//  * n: power constant
// =============================================================================
inline double R_vm_martynov(const double Q2, const double Mv, const double c,
                            const double n) {
  const double Mv2 = Mv * Mv;
  return pow((c * Mv2 + Q2) / (c * Mv2), n) - 1.;
}

// =============================================================================
// Dipole form factor to relate real photo-production cross section to sigma_T
//
// This form deviates from the classical VMD form (which has a fixed power of
// 2), to better fit the world data for rho0 production.
//
// Cf.:
//      A. Airapetian et al, "Exclusive Leptoproduction of rho0 Mesons on
//      Hydrogen at Intermediate W Values", EPJ C 17 (2000) 389-398
//
//      Adams et al., "Diffractive production of ρ0 mesons in muon–proton
//      interactions 470 GeV", ZPC74 (1997) 237-261.
//
// An optimized power of n=2.575 was taken from a fit result from
//
//      M Tytgat, "Diffractive production of ρ0 and ω vector mesons at HERMES",
//      DESY-Thesis 2001-018 (2001)
//
// Arguments:
//  * Q2: -photon mass squared
//  * Mv: VM mass
//  * n: power constant (2 for VMD; 2.575 from HERMES fit)
// =============================================================================
inline double dipole_ff_vm(const double Q2, const double Mv, const double n) {
  const double Mv2 = Mv * Mv;
  return pow(Mv2 / (Mv2 + Q2), n);
}

// =============================================================================
// t-channel photo-production cross section for heavy vector mesons
// Values are in units of nb/GeV^2.
//
// Formulism from
//      brodsky et. al., Phys.Lett.B498:23-28,2001
//      (http://arxiv.org/abs/hep-ph/0010343)
//
// Both dsigma/dt and dsigma/d(exp_bt) are available. The latter is useful for a
// more efficient MC implementation
//
// Parameter values from Eric's fit to the J/psi world data:
//  * c2g: 6.499e3 [GeV^-2]
//  * c3g: 2.894e3 [GeV^-2]
//  * b: 1.13 [GeV^-2]
//
// Arguments:
//  * s: photon-target system invariant mass squared
//  * t: mandelstam t
//  * Mt: target mass
//  * Mv: VM mass
//  * b: b-parameter of target form-factor
//  * c2g: 2-gluon constant
//  * c3g: 3-gluon constant
//
// =============================================================================
inline double dsigma_dt_vm_brodsky(const double s, const double t,
                                   const double Mt, const double Mv,
                                   const double b, const double c2g,
                                   const double c3g = 0) {
  const double Mt2 = Mt * Mt;
  const double Mv2 = Mv * Mv;
  const double x = (2. * Mt * Mv + Mv2) / (s - Mt2);
  // form factor
  const double ff = exp(b * t);
  // phase space factor
  const double v = 1. / (16. * TMath::Pi());
  // 2 gluon term
  const double A2g = c2g * v * (1 - x) * (1 - x) / Mv2;
  // 3 gluon term
  const double A3g = c3g * v / (Mv2 * Mv2);
  return (A2g + A3g) * ff;
}
inline double dsigma_dexp_bt_vm_brodsky(const double s, const double Mt,
                                        const double Mv, const double b,
                                        const double c2g,
                                        const double c3g = 0) {
  const double Mt2 = Mt * Mt;
  const double Mv2 = Mv * Mv;
  const double x = (2. * Mt * Mv + Mv2) / (s - Mt2);
  // form factor (absorbed by the jacobian for d(bt) -> d(exp_bt)
  const double ff = 1.;
  // phase space factor
  const double v = 1. / (16. * TMath::Pi());
  // 2 gluon term
  const double A2g = c2g * v * (1 - x) * (1 - x) / Mv2;
  // 3 gluon term
  const double A3g = c3g * v / (Mv2 * Mv2);
  // extra jacobian for dt -> d(bt)
  const double jacobian = 1 / b;
  return (A2g + A3g) * ff * jacobian;
}

} // physics
} // lager

#endif
