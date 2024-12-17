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

#ifndef LAGER_PHYSICS_VM_LOADED
#define LAGER_PHYSICS_VM_LOADED

#include <TMath.h>
#include <cmath>

#include <lager/physics/kinematics.hh>

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
// Dipole (multipole) form factor to relate real photo-production cross section
// to sigma_T
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
inline double multipole_ff_vm(const double Q2, const double Mv,
                              const double n) {
  return dipole_ff_vm(Q2, Mv, n);
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

// =============================================================================
// Electroproduction of phi mesons according the formalism
// of the CLAS12 exclusive phi meson electroproduction proposal
// https://www.jlab.org/exp_prog/proposals/12/PR12-12-007.pdf
//
// This is the transverse part of the cross section, also valid for
// photoproduction. Note that this is the integrated cross section, which needs
// to be multiplied with a normalized form factor F(t)/F_int for a differential
// cross section
// =============================================================================
inline double sigmaT_phi_clas(const double Q2, const double W, const double Mt,
                              const double Mv, const double alpha_1,
                              const double alpha_2, const double alpha_3,
                              const double nu_T) { //, const double B0,
  const double Wth = Mt + Mv;
  const double cT =
      alpha_1 * pow((1 - Wth * Wth / (W * W)), alpha_2) * pow(W, alpha_3);
  const double sigmaT = cT * multipole_ff_vm(Q2, Mv, nu_T);
  return sigmaT;
}
inline double R_phi_clas(const double Q2, const double Mv, const double c_R) {
  const double Mv2 = Mv * Mv;
  return c_R * Q2 / Mv2;
}

  
// =============================================================================
// Electroproduction of phi mesons according the GPD formalism
// of Hatta, Passek, Schoenleber
//
// This is the transverse part of the cross section, also valid for
// photoproduction. Note that this is the integrated cross section, which needs
// to be multiplied with a normalized form factor F(t)/F_int for a differential
// cross section
// =============================================================================
inline double sigmaT_phi_hatta(const double Q2, const double W, const double Mt,
                              const double Mv, const double alpha_1,
                              const double alpha_2, const double alpha_3,
                              const double nu_T) { //, const double B0,                                                                                                                                                                                                                  
  const double Wth = Mt + Mv;
  double cT_temp = 0.0;
  double cT_clas12 = alpha_1 * pow((1 - Wth * Wth / (W * W)), alpha_2) * pow(W, alpha_3);
  double cT_hatta = (140.0/pow(sqrt(Q2), 9.4)) * pow(W*W-Mv*Mv, 0.71) * pow(W,3.8); //Nominally, a1 = 0.045, a2 = 0.866, a3 = 4.177                                                                                                                                                      
    //alpha_1 * pow((2/sqrt(Q2)), alpha_2) * pow(W, alpha_3); //Nominally, a1 = 0.045, a2 = 0.866, a3 = 4.177                                                                                                                                                                            
  //std::cout << " Q2 = " << Q2 << "Wth = " << Wth <<  " W = " << W << std::endl;
  //std::cout << "0.045 * pow((2/sqrt(Q2)), 9.65) = " << 0.045 * pow((2/sqrt(Q2)), 9.65) << " pow(W*W-Mv*Mv, 0.866) = " << pow(W*W-Mv*Mv, 0.866) << " pow(W,4.177) = " << pow(W,4.177) << std::endl;
  //std::cout << "CLAS ct = " << cT_clas12 << " hatta cT = " << cT_hatta << " ratio = " << cT_clas12/cT_hatta << std::endl;

  //Choose whichever is smaller between Hatta and CLAS12 models                                                                                                                                                                                                                          
  if(cT_clas12 > cT_hatta) cT_temp = cT_hatta;
  else cT_temp = cT_clas12;

  const double cT = cT_temp;
  //std::cout << "final cT = " << cT << std::endl;
  const double sigmaT = cT * multipole_ff_vm(Q2, Mv, nu_T);
  return sigmaT;
}

  inline double R_phi_hatta(const double Q2, const double Mv, const double c_R) {
  const double Mv2 = Mv * Mv;
  return c_R * Q2 / Mv2;
}

inline double exp_ff_normalized(const double Q2, const double W, const double t,
                                const double Mt, const double Mv,
                                const double B0, const double alphaP) {
  const double B = B0 + 4 * alphaP * std::log(W);
  const double F = exp(B * t);
  const double t_min = t_range(W * W, Q2, Mt, Mv, Mt).max;
  const double F_int = exp(B * t_min) / B;
  return F / F_int;
}
inline double dipole_ff_normalized(const double Q2, const double W,
                                   const double t, const double Mt,
                                   const double Mv, const double Mg2) {
  const double Mg8 = pow(Mg2, 4);
  const double F = Mg8 / pow(Mg2 - t, 4);
  const double t_min = t_range(W * W, Q2, Mt, Mv, Mt).max;
  const double F_int = Mg8 / (3 * pow(Mg2 - t_min, 3));
  return F / F_int;
}
// =============================================================================
// J/psi production in the holographic model by Mamo & Zahed
//
// This is the exact formalism used in https://inspirehep.net/literature/2110821
// =============================================================================
inline double F_kinematic_holographic(double s, double t, double Mt,
                                      double Mv) {
  double F =
      (1 / (4096 * Mv * Mv)) *
      (-9 * pow(Mv, 10) + pow(Mv, 8) * (-32 + 68 * Mt * Mt + 28 * s + 37 * t) +
       2 * pow(Mv, 6) *
           (256 * pow(Mt, 4) + 8 * Mt * Mt * (32 * s - 3 * t) +
            t * (56 - 40 * s - 29 * t)) +
       2 * pow(Mv, 4) *
           (-136 * pow(Mt, 6) + 64 * s * s - 56 * pow(s, 3) +
            8 * pow(Mt, 4) * (8 + 27 * s - 64 * t) + 3 * t * t * (-24 + 7 * t) +
            4 * s * t * (-4 + 9 * t) -
            4 * Mt * Mt *
                (6 * s * s + 32 * s * (1 + 4 * t) + t * (-4 + 25 * t))) +
       Mv * Mv *
           (144 * pow(Mt, 8) + 144 * pow(s, 4) - 192 * s * s * t +
            96 * pow(s, 3) * t - 16 * s * (-4 + t) * t * t +
            (80 - 13 * t) * pow(t, 3) + 96 * pow(Mt, 6) * (-6 * s + 7 * t) +
            32 * pow(Mt, 4) * (27 * s * s - 6 * t - 39 * s * t + 8 * t * t) +
            16 * Mt * Mt *
                (-36 * pow(s, 3) + 30 * s * s * t + 24 * s * t * (1 + 2 * t) +
                 t * t * (-4 + 17 * t))) -
       t * (2 * Mt * Mt - 2 * s - t) *
           (64 * pow(Mt, 4) + 8 * pow(Mt, 6) - 8 * pow(s, 3) +
            76 * pow(Mt, 4) * t - 16 * t * t - 90 * Mt * Mt * t * t +
            pow(t, 3) + 4 * s * s * (16 + 6 * Mt * Mt + 3 * t) -
            2 * s * (12 * pow(Mt, 4) + 3 * t * t + Mt * Mt * (64 + 44 * t))));
  return F;
}
// Arguments:
//  - Q2, W, t  --> kinematics
//  - Mt, Mv    --> target mass and VM mass
//  - A0, m_A   --> tripole parameters for the A GFF
//  - C0, m_C   --> tripole parameters for the C GFF
//  - N         --> Additional free normalization constant
inline double dsigma_dt_holographic(const double Q2, const double W,
                                    const double t, const double Mt,
                                    const double Mv, const double A0,
                                    const double m_A, const double C0,
                                    const double m_C, const double N) {
  const double s = W * W;
  const double abst = -t;
  const double q_gamma =
      (1 / (2 * pow(s, 0.5))) *
      pow(s * s - 2 * (-Q2 + Mt * Mt) * s + (-Q2 - Mt * Mt) * (-Q2 - Mt * Mt),
          0.5);
  const double q_jpsi = (1 / (2 * pow(s, 0.5))) *
                        pow(s * s - 2 * (Mv * Mv + Mt * Mt) * s +
                                (Mv * Mv - Mt * Mt) * (Mv * Mv - Mt * Mt),
                            0.5);
  const double E_jpsi = pow(Mv * Mv + q_jpsi * q_jpsi, 0.5);
  const double E_gamma = pow(-Q2 + q_gamma * q_gamma, 0.5);
  const double const_h =
      1 / (64 * Mt * Mt * pow(s - Mt * Mt, 2) * 4 * TMath::Pi());
  const double eta = (Mv * Mv) / (2 * (s - Mt * Mt) - Mv * Mv - abst);
  const double Ak = A0 / pow(1 + abst / (m_A * m_A), 3);
  const double Ck = C0 / pow(1 + abst / (m_C * m_C), 3);
  const double F = F_kinematic_holographic(s, -abst, Mt, Mv);
  const double dsigma_dt =
      N * N * const_h *
      (Ak * Ak + 2 * eta * eta * Ak * 4 * Ck + pow(eta, 4) * pow(4 * Ck, 2)) /
      (A0 * A0) * F * (2 * abst + 8 * Mt * Mt);
  return dsigma_dt;
}

// =============================================================================
// General VM decay distributions in the VM helicity frame for the
// cases of  (1) VM --> Scaler+scaler
//       and (2) VM -> fermion+fermion
// Note that case (1) corresponds equation 31 (for W0) in
//     K. Schilling et al, Nucl.Phys.B 15 (1970) 397-412
// while case (2) has the corresponding formula for decay in spin-1/2
// particles
//
// Both expressions are a function of the decay angles in the VM helicity
// frame
//  cth: cosine of the polar angle
//  phi: azimuthal angle with the VM production plane
//  sdme_04_00:  r^04_00 (=rho^0_00 for photoproduction)
//  sdme_04_10:  Re(r^04_10) (=Re(rho^0_10) for photoproduction)
//  sdme_04_1m1: r^04_1-1 (=rho^0_1-1 for photoproduction)
// =============================================================================
inline double vm_decay_scalars(const double cth, const double phi,
                               const double sdme_04_00, const double sdme_04_10,
                               const double sdme_04_1m1) {
  const double theta = acos(cth);
  const double sth = sin(theta);
  const double factor = 3 / (4. * TMath::Pi());
  const double t1 = 0.5 * (1 - sdme_04_00);
  const double t2 = 0.5 * (3 * sdme_04_00 - 1) * cth * cth;
  const double t3 =
      -1 * TMath::Sqrt2() * sdme_04_10 * sin(2 * theta) * cos(phi);
  const double t4 = -1 * sdme_04_1m1 * sth * sth * cos(2 * phi);
  return t1 + t2 + t3 + t4;
}
inline double vm_decay_fermions(const double cth, const double phi,
                                const double sdme_04_00,
                                const double sdme_04_10,
                                const double sdme_04_1m1) {
  const double theta = acos(cth);
  const double sth = sin(theta);
  const double factor = 3 / (4. * TMath::Pi());
  const double t1 = 0.5 * (1 + sdme_04_00);
  const double t2 = -0.5 * (3 * sdme_04_00 - 1) * cth * cth;
  const double t3 = 1 * TMath::Sqrt2() * sdme_04_10 * sin(2 * theta) * cos(phi);
  const double t4 = 1 * sdme_04_1m1 * sth * sth * cos(2 * phi);
  LOG_JUNK2("vm_decay_fermions",
            "cos(theta): " + std::to_string(cth) +
                ", theta: " + std::to_string(theta) +
                ", t1: " + std::to_string(t1) + ", t2: " + std::to_string(t2) +
                ", t3: " + std::to_string(t3) + ", t4: " + std::to_string(t4));
  return t1 + t2 + t3 + t4;
}

} // namespace physics
} // namespace lager

#endif
