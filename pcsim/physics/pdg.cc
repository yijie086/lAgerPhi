#include <pcsim/physics/pdg.hh>

namespace pcsim {
namespace physics {

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
