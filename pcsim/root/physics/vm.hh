#ifndef PCSIM_ROOT_PHYSICS_VM_LOADED
#define PCSIM_ROOT_PHYSICS_VM_LOADED

// ROOT module with TF1 versions of various physics functions

#include <Math/Vector4D.h>
#include <TF1.h>
#include <pcsim/core/interval.hh>

namespace pcsim {
namespace root {
namespace physics {
namespace vm {

using range_type = interval<double>;

// epsilon as a function of Q2 (W and and beam energy are parameters)
//
// TF1 parameters *x and *par
// x[0]: Q2
// par[0]: W
// par[1]: Ebeam
TF1* epsilon_Q2(const range_type& Q2lim);

// =============================================================================
// t-channel photo-production cross section for heavy vector mesons
// Values are in units of nb/GeV^2.
//
// Formulism from
//      brodsky et. al., Phys.Lett.B498:23-28,2001
//      (http://arxiv.org/abs/hep-ph/0010343)
//
// =============================================================================
//
// TF1 parameters *x and *par
// x[0]: t (sign force-added, in case we want to draw -t)
// par[0]: 2-gluon amplitude
// par[1]: 3-gluon amplitude
// par[2]: b
// par[3]: VM mass
// par[4]: W
TF1* dsigma_dt_vm_brodsky(const range_type& tlim);

// same but differential in d/d(exp(bt)) for more efficient integration
//
// TF1 parameters *x and *par
// x[0]: t (unused)
// par[0]: 2-gluon amplitude
// par[1]: 3-gluon amplitude
// par[2]: b
// par[3]: VM mass
// par[4]: W
TF1* dsigma_dexp_bt_vm_brodsky(const range_type& tlim);

// t-integrated cross section as a function of W
//
// TF1 parameters *x and *par
// x[0]: W
// par[0]: 2-gluon Amplitude
// par[1]: 3-gluon Amplitude
// par[2]: b
// par[3]: VM mass
TF1* sigma_vm_brodsky_W(const range_type& Wlim);

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
// =============================================================================
//
// TF1 parameters *x and *par
// x[0]: Q2
// par[0]: Mv
// par[1]: c
// par[2]: n
TF1* R_vm_martynov(const range_type& Q2lim);

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
// =============================================================================
//
// TF1 parameters *x and *par
// x[0]: Q2
// par[0]: Mv
// par[1]: n
TF1* dipole_ff_vm(const range_type& Q2lim);

// =============================================================================
// t-integrated VM cross section as a function of (Q2, W) (x-val: W)
// =============================================================================
//
// TF1 parameters *x and *par
// x[0]: W
// par[0]: 2-gluon Amplitude
// par[1]: 3-gluon Amplitude
// par[2]: b
// par[3]: VM mass
// par[4]: Q2
// par[5]: W (unused)
// par[6]: Ebeam
// par[7]: dipole_n
// par[8]: R_c
// par[9]: R_n
//
// TF1 parameters *x and *par
TF1* sigma_vm_brodsky_W_Q2(const range_type& Wlim);

// =============================================================================
// VM cross section differential in t as a function of (Q2, W) (x-val: t)
// =============================================================================
//
// TF1 parameters *x and *par
// x[0]: t
// par[0]: 2-gluon Amplitude
// par[1]: 3-gluon Amplitude
// par[2]: b
// par[3]: VM mass
// par[4]: Q2
// par[5]: W
// par[6]: Ebeam
// par[7]: dipole_n
// par[8]: R_c
// par[9]: R_n
TF1* dsigma_dt_vm_brodsky_Q2W(const range_type& tlim);

// =============================================================================
// r00_04 matrix element
// =============================================================================
//
// TF1 parameters *x and *par
// x[0]: Q2
// par[0]: W
// par[1]: Mv
// par[2]: Ebeam
// par[3]: R_c
// par[4]: R_n
TF1* r00_04(const range_type& Q2lim);

// =============================================================================
// SCHC theta angle (cos(theta))
// =============================================================================
//
// TF1 parameters *x and *par
// x[0]: cos(theta)
// par[0]: Q2
// par[1]: W
// par[2]: Mv
// par[3]: Ebeam
// par[4]: R_c
// par[5]: R_n
TF1* ctheta_schc(const range_type& cthlim = {-1, 1});

// =============================================================================
// t-integrated VM cross section differential in ctheta as a function of (Q2, W)
// =============================================================================
//
// TF1 parameters *x and *par
// x[0]: cos(theta)
// par[0]: 2-gluon Amplitude
// par[1]: 3-gluon Amplitude
// par[2]: b
// par[3]: VM mass
// par[4]: Q2
// par[5]: W
// par[6]: Ebeam
// par[7]: dipole_n
// par[8]: R_c
// par[9]: R_n
TF1* dsigma_dctheta_vm_brodsky_Q2W(const range_type& cthlim = {-1, 1});

} // namespace vm
} // namespace physics
} // namespace root
} // namespace pcsim

#endif
