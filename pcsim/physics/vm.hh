#ifndef PCSIM_PHYSICS_VM_LOADED
#define PCSIM_PHYSICS_VM_LOADED

#include <TMath.h>
#include <cmath>
#include <pcsim/core/particle.hh>

// =============================================================================
// VM Physics routines
// =============================================================================

namespace pcsim {
namespace physics {

// =============================================================================
// R VM parameterization 
//
// from 
//      Martynov, et. al., “Photoproduction of Vector Mesons in the Soft Dipole
//      Pomeron Model.” PRD 67 (7), 2003. doi:10.1103/PhysRevD.67.074023.
// eq. 31.
//
// J/psi parameters a and n 
//  * a: 2.164
//  * n: 2.131
// from eq. 18:
//
//      R. Fiore et al., "Exclusive Jpsi electroproduction in a dual model."
//      PRD80:116001, 2009"
// =============================================================================
inline double R_vm_martynov(const double Q2, const particle& vm,
                            const double a, const double n) {
  const double Mv2 = vm.mass * vm.mass;
  return pow((a * Mv2 + Q2) / (a * Mv2), n) - 1.;
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
// =============================================================================
inline double dipole_ff_vm_hermes(const double Q2, const particle& vm,
                                  const double n) {
  const double Mv2 = vm.mass * vm.mass;
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
// =============================================================================
inline double dsigma_dt_vm_brodksy(const double s, const double t,
                                   const particle& target, const particle& vm,
                                   const double b, const double c2g,
                                   const double c3g = 0) {
  const double Mt = target.mass;
  const double Mt2 = Mt * Mt;
  const double Mv = vm.mass;
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
inline double dsigma_dexp_bt_brodksy(const double s, const double t,
                                     const particle& target, const particle& vm,
                                     const double b, const double c2g,
                                     const double c3g = 0) {
  const double Mt = target.mass;
  const double Mt2 = Mt * Mt;
  const double Mv = vm.mass;
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
} // pcsim

#endif
