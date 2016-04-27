#ifndef PCSIM_PHYSICS_PDG_LOADED
#define PCSIM_PHYSICS_PDG_LOADED

#include <cstdint>
#include <memory>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

namespace pcsim {
namespace physics {

constexpr int32_t pdg_nuclear_id(unsigned A, unsigned Z) {
  return 99000000 + A * 1000 + Z;
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
  J_Psi = 443,
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
  // unknown
  unknown = -9999
};

// wrapper around ROOTs PDG database that automatically adds info on common
// nuclei (therefor avoiding the need for a custom database)
// implemented using a single persistent static database
// note: nuclear masses in GeV from
// http://hyperphysics.phy-astr.gsu.edu/hbase/pertab
inline TParticlePDG* pdg_particle(const int32_t id) {
  // use a leaky naked pointer to the TDatabasePDG, as using a shared_ptr
  // segfaults on destruction
  static TDatabasePDG* db;
  if (!db) {
    db = new TDatabasePDG{};
    db->AddParticle("H-2", "Deuteron", 1.875613, true, 0, 1, "nucleus",
                    static_cast<int32_t>(pdg_id::H2));
    db->AddParticle("H-3", "Triton", 2.808921, true, 0, 1, "nucleus",
                    static_cast<int32_t>(pdg_id::H3));
    db->AddParticle("He-3", "Helium-3", 2.808391, true, 0, 1, "nucleus",
                    static_cast<int32_t>(pdg_id::He3));
    db->AddParticle("He-4", "Helium-4", 3.727379, true, 0, 1, "nucleus",
                    static_cast<int32_t>(pdg_id::He4));
    db->AddParticle("C-12", "Carbon-12", 11.1750, true, 0, 1, "nucleus",
                    static_cast<int32_t>(pdg_id::C12));
    db->AddParticle("N-14", "Nitrogen-14", 13.0403, true, 0, 1, "nucleus",
                    static_cast<int32_t>(pdg_id::N14));
  }
  return db->GetParticle(id);
}
inline TParticlePDG* pdg_particle(const pdg_id id) {
  return pdg_particle(static_cast<int32_t>(id));
}

// some often used particles (add as needed, don't forget pdg.cc)
extern const TParticlePDG& PDG_ELECTRON;
extern const TParticlePDG& PDG_PROTON;
extern const TParticlePDG& PDG_ANTIPROTON;
extern const TParticlePDG& PDG_PI_PLUS;
extern const TParticlePDG& PDG_PI_MINUS;
extern const TParticlePDG& PDG_K_PLUS;
extern const TParticlePDG& PDG_K_MINUS;
extern const TParticlePDG& PDG_PHOTON;
extern const TParticlePDG& PDG_JPSI;

} // physics
} // pcsim

#endif
