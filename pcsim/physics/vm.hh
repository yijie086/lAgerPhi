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
double R_vm_martynov(const double Q2, const particle& vm, const double R_a,
                     const double R_n) {
  const double Mv2 = vm.mass * vm.mass;
  return pow((R_a * Mv2 + Q2) / (R_a * Mv2), R_n) - 1.;
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
// An optimized power of 2.575 was taken from a fit result from
//
//      M Tytgat, "Diffractive production of ρ0 and ω vector mesons at HERMES",
//      DESY-Thesis 2001-018 (2001)
// =============================================================================
double dipole_ff_vm(const double Q2, const particle& vm,
                    const double dipole_n) {
  const double Mv2 = vm.mass * vm.mass;
  return pow(Mv2 / (Mv2 + Q2), dipole_n);
}

} // physics
} // pcsim

#endif
