#include <pcsim/physics/pdg.hh>

namespace pcsim {
namespace physics {

TParticlePDG* pdg_particle(const int32_t id) {
  // implemented using a single persistent static database
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

// some often used particles (add as needed, don't forget the header)
const TParticlePDG& PDG_ELECTRON{*pdg_particle(pdg_id::e_minus)};
const TParticlePDG& PDG_PROTON{*pdg_particle(pdg_id::p)};
const TParticlePDG& PDG_ANTIPROTON{*pdg_particle(pdg_id::p_bar)};
const TParticlePDG& PDG_PI_PLUS{*pdg_particle(pdg_id::pi_plus)};
const TParticlePDG& PDG_PI_MINUS{*pdg_particle(pdg_id::pi_minus)};
const TParticlePDG& PDG_K_PLUS{*pdg_particle(pdg_id::K_plus)};
const TParticlePDG& PDG_K_MINUS{*pdg_particle(pdg_id::K_minus)};
const TParticlePDG& PDG_PHOTON{*pdg_particle(pdg_id::gamma)};
const TParticlePDG& PDG_JPSI{*pdg_particle(pdg_id::J_Psi)};

} // physics
} // pcsim
