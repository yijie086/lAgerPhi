#ifndef PCSIM_CORE_PDG_LOADED
#define PCSIM_CORE_PDG_LOADED

#include <TDecayChannel.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <cstdint>
#include <memory>

namespace pcsim {

// Pentaquark particle codes
// (slight deviation from PDG scheme:
//   - parity code to allow for different parities in case exact parity isn't
//     known
//   - extra code in case we don't know the exact mass assignments of the
//     penaquarks

#if __cplusplus < 201402L
// extra utility function needed pre-C++14
// in C++14 and later we can just calculate this in the constexpr function and
// store it as a variable
constexpr int pdg_pentaquark_id_pcode(int parity) {
  return (parity == -1) ? 2 : ((parity == 1) ? 1 : 0);
}
#endif

constexpr int32_t pdg_pentaquark_id(int q1, int q2, int q3, int q4, int qbar,
                                    int twoJ, int parity = 0,
                                    int extracode = 0) {
#if __cplusplus < 201402L
  return extracode * 100000000 + pdg_pentaquark_id_pcode(parity) * 10000000 +
         9 * 1000000 + q1 * 100000 + q2 * 10000 + q3 * 1000 + q4 * 100 +
         qbar * 10 + twoJ + 1;
#else
  int pcode = 0;
  if (parity == -1) {
    pcode = 2;
  } else if (parity == 1) {
    pcode = 1;
  }
  return extracode * 100000000 + pcode * 10000000 + 9 * 1000000 + q1 * 100000 +
         q2 * 10000 + q3 * 1000 + q4 * 100 + qbar * 10 + twoJ + 1;
#endif
}

// add possibility for nuclear id
// prefixed by 4 9's
constexpr int32_t pdg_nuclear_id(unsigned A, unsigned Z) {
  return 9999000000 + A * 1000 + Z;
}

// particle ID enum
enum class pdg_id : int32_t {
  // partons
  d = 1,
  d_bar = -1,
  u = 2,
  u_bar = -2,
  s = 3,
  s_bar = -3,
  c = 4,
  c_bar = -4,
  b = 5,
  b_bar = -5,
  t = 6,
  t_bar = -6,
  // leptons
  e_minus = 11,
  e_plus = -11,
  nu_e = 12,
  nu_e_bar = -12,
  mu_minus = 13,
  mu_plus = -13,
  nu_mu = 14,
  nu_mu_bar = -14,
  tau_minus = 15,
  tau_plus = -15,
  nu_tau = 16,
  nu_tau_bar = -16,
  // bosons
  g = 21,
  gamma = 22,
  Z_0 = 23,
  W_plus = 24,
  W_minus = -24,
  H_0 = 25,
  reggeon = 28,
  pomeron = 29,
  // scalar mesons
  pi_0 = 111,
  pi_plus = 211,
  pi_minus = -211,
  eta = 221,
  K_0 = 311,
  K_0_bar = -311,
  K_L_0 = 130,
  K_S_0 = 310,
  K_plus = 321,
  K_minus = -321,
  eta_prime = 331,
  // vector mesons
  rho_0 = 113,
  rho_plus = 213,
  rho_minus = -213,
  omega = 223,
  K_star_0 = 313,
  K_star_0_bar = -313,
  K_star_plus = 323,
  K_star_minus = -323,
  phi = 333,
  J_psi = 443,
  upsilon = 553,
  // baryons
  n = 2112,
  n_bar = -2122,
  p = 2212,
  p_bar = -2212,
  // nuclei
  H1 = 2212, // same as proton
  H2 = pdg_nuclear_id(2, 1),
  H3 = pdg_nuclear_id(3, 1),
  He3 = pdg_nuclear_id(3, 2),
  He4 = pdg_nuclear_id(4, 2),
  C12 = pdg_nuclear_id(12, 6),
  N14 = pdg_nuclear_id(14, 7),
  // LHCb pentaquark hypetheses (in order of likelihood for the 4450 and 4380)
  //                              q1 q2 q3 q4 qb 2J P  Extra
  Pc_4450_52p = pdg_pentaquark_id(4, 2, 2, 1, 4, 5, 1, 1),
  Pc_4450_52m = pdg_pentaquark_id(4, 2, 2, 1, 4, 5, -1, 2),
  Pc_4450_32p = pdg_pentaquark_id(4, 2, 2, 1, 4, 3, 1, 3),
  Pc_4450_32m = pdg_pentaquark_id(4, 2, 2, 1, 4, 3, -1, 4),
  Pc_4380_32m = pdg_pentaquark_id(4, 2, 2, 1, 4, 3, -1, 1),
  Pc_4380_32p = pdg_pentaquark_id(4, 2, 2, 1, 4, 3, 1, 2),
  Pc_4380_52m = pdg_pentaquark_id(4, 2, 2, 1, 4, 5, -1, 3),
  Pc_4380_52p = pdg_pentaquark_id(4, 2, 2, 1, 4, 5, 1, 4),
  // unknown
  unknown = -9999
};

// Get PID info from the buildin ROOT PDG database
// note: the custom particles (nuclei, pentaquarks, ...) under pdg_id are added
// to the database.
// Feel free to add other particles as-needed to the database in pdg.cc
TParticlePDG* pdg_particle(const pdg_id id);

} // ns pcsim

#endif
