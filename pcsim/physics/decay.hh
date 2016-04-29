#ifndef PCSIM_PHYSICS_DECAY_LOADED
#define PCSIM_PHYSICS_DECAY_LOADED

#include <TLorentzVector.h>
#include <TRandom.h>

namespace pcsim {
namespace physics {

// J/Psi leptonic decay, assuming helicity conservation (1+cos^2\theta_CM).
// arguments:
//  rng: the random generator to be used
//  jpsi: J/Psi 4-vector
//  ml: lepton mass (electron or muon mass)
//  lplus, lminus: reference to the decay lepton 4-vectors
void decay_jpsi_lepton(std::shared_ptr<TRandom> rng, const TLorentzVector& jpsi,
                       const double ml, TLorentzVector& lplus,
                       TLorentzVector& lminus);
// Pc decay into p+J/Psi, ignoring the J/Psi spin and parity (using an isotropic
// decay distribution)
void decay_pc_iso(std::shared_ptr<TRandom> rng, const TLorentzVector& pc,
                  TLorentzVector& proton, TLorentzVector& jpsi);
} // physics
} // pcsim

#endif
