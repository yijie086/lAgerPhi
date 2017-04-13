#ifndef PCSIM_CORE_PARTICLE_LOADED
#define PCSIM_CORE_PARTICLE_LOADED

#include <Math/LorentzRotation.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

#include <TRandom.h>
#include <memory>
#include <pcsim/core/interval.hh>
#include <pcsim/core/pdg.hh>

namespace pcsim {

// =============================================================================
// particle 
//
// cornerstone class containing particle info
// =============================================================================
class particle {
public:
  using ROOT::Math::XYZTVector::BetaVector;
  using ROOT::Math::XYZTVector;
  using ROOT::Math::Boost;
  using ROOT::Math::XYZVector;
  using ROOT::Math::Polar3DVector;

  // particle status
  enum class status_code : int {
    BEAM = 11,
    SECONDARY_BEAM = 12,
    SCAT = 13,
    UNSTABLE = 21,
    UNSTABLE_SCHC = 22,
    FINAL = 30,
    DECAYED = 31,
    DECAYED_SCHC = 32,
    OTHER = 99
  };

  // constructors
  particle() = default;
  particle(const particle& rhs) = default;
  // particle with a given id
  particle(const pdg_id id);
  // particle with a given momentum 3-vector
  particle(const pdg_id id, const XYZVector& p3);
  // particle in a given direction with a given energy
  particle(const pdg_id id, const XYZVector& direction, const double E);
  // particle with a given 4-momentum (can be off-shell)
  particle(const pdg_id id, const XYZTVector& p);
  // RNG constructors usefull in particular for unstable particles with a mass
  // and lifetime that needs to be generated
  particle(const pdg_id id, std::shared_ptr<TRandom> rng);
  particle(const pdg_id id, const XYZVector& p3, std::shared_ptr<TRandom> rng);

  // particle ID code
  pdg_id type() const { return type_; }
  int type() const { return static_cast<int>(type_); }

  // full PDG info from the TDatabasePDG
  TParticlePDG* pdg() { return pdg_; }

  // particle status
  status_code status() const { return status_; }
  int status() const { return static_cast<int>(status_); }
  void update_status(const status_code ns) { status_ = ns; }

  // particle properties
  int charge() const { return charge_; }
  double mass() const { return mass_; }
  double mass2() const { return mass_ * mass_; }
  double width() const { return width_; }
  double lifetime() const { return lifetime_; }
  interval<int> parent() const { return parent_; }
  interval<int> daughter() const { return daughter_; }
  int daughter_begin() const { return daughter_.min; }
  int daughter_end() const { return daughter_.max; }
  int n_daughters() const { return daughter_.max - daughter_.min; }
  int parent_first() const { return parent_.min; }
  int parent_second() const { return parent_.max; }
  int n_parents() const { return parent_.max - parent_.min; }
  double momentum() const { return p_.Vect().Mag(); }
  double energy() const { return p_.E(); }

  // get a reference to the momentum/vertex 4-vector
  XYZTVector& p() { return p_; };
  const XYZTVector& p() const { return p_; };
  XYZTVector& vertex() { return vertex_; }
  const XYZTVector& vertex() const { return vertex_; }

  // transformations etc.
  void boost(const BetaVector& bv);
  void boost(const Boost& b) { p = b * p; }
  // rotate from a coordate system where v moved along the z-axis, to the
  // coordate system of v
  // (algorithm taken from TVector3::RotateUz)
  void rotate_uz(const XYZVector& v);
  void rotate_uz(const XYZTVector& v);
  void rotate_uz(const particle& pv);

  // add a daughter track
  void add_daughter(const int index);
  // add a parent track
  void add_parent(const int index);
  // set the parent tracks
  void set_parents(const interval<int> indices) { parent_ = indices; }

private:
  void set_mass_lifetime(std::shared_ptr<TRandom> rng);
  XYZTVector make_4vector(const XYZVector v, const double t) {
    return {v.X(), v.Y(), v.Z(), v.t()};
  }

  pdg_id type_{pdg_id::UNKNOWN};
  TParticlePDG* pdg_{nullptr};
  int charge_{0};
  double width_{0};
  status_code status_{status_code::OTHER};

  // actual mass and lifetime for this particle, can deviate from pole values
  // for unstable particles!
  double mass_{0};
  double lifetime_{0};
  XYZTVector p_{0, 0, 0, 0};
  XYZTVector vertex_{0, 0, 0, 0};
  // parent indices store the first (and optional second) parent of a
  // particle. -1 if not stored
  interval<int> parent_{-1, -1};
  // daughter indices are encoded from [begin, end) where begin is the
  // first index and end one past the last index. (-1, -1) if none
  interval<int> daughter_{-1, -1};
};

} // namespace pcsim

// =============================================================================
// Particle implementation
//
// This is a work-horse class called, and therefore all its members are defined
// inline
// =============================================================================
namespace pcsim {

//
// CONSTRUCTORS
//
// particle with a given id
inline particle::particle(const pdg_id)
    : type_{id}
    , pdg_{pdg_particle(id)}
    , charge_{static_cast<int>(pdg_->Charge() / 3)}
    , width_{pdg_->Width()}
    , status_{width_ > 0 ? status_code::UNSTABLE : status_code::FINAL} {}
// particle with a given momentum 3-vector
inline particle::particle(const pdg_id id, const XYZVector& p3)
    : particle{id}
    , mass_{pdg->Mass()}
    , p_{make_4vector(p3, sqrt(mass_ * mass_ + p_.Mag2()))} {}
// particle in a given direction with a given energy
inline particle::particle(const pdg_id id, const XYZVector& direction,
                          const double E)
    : particle{id}
    , mass_{pdg->Mass()}
    , p_{make_4vector(direction.Unit() * sqrt(E * E - mass_ * mass_), E)} {}
// particle with a given 4-momentum (can be off-shell)
inline particle::particle(const pdg_id id, const XYZTVector& p)
    : particle{id}, mass_{p.M()}, p_{p} {}
// RNG constructors usefull in particular for unstable particles with a mass
// and lifetime that needs to be generated
inline particle::particle(const pdg_id id, std::shared_ptr<TRandom> rng)
    : particle(id) {
  set_mass_lifetime(std::move(rng));
}
inline particle::particle(const pdg_id id, const XYZVector& p3,
                          std::shared_ptr<TRandom> rng)
    : particle{id} {
  set_mass_lifetime(std::move(rng));
  p_.SetXYZT(p.X(), p.Y(), p.Z(), sqrt(mass_ * mass_ + p.Mag2()));
}
//
// parent/daughter modifiers
//
inline void particle::add_daughter(const int index) {
  if (daughter_.min > index || daughter_.min < 0) {
    daughter_.min = index;
  }
  if (daughter_.max >= index || daughter_.max < 0) {
    daughter_.max = index + 1;
  }
}
inline void particle::add_parent(const int index) {
  if (parent_.min < 0) {
    parent_.min = index;
  } else if (parent_.max < 0) {
    parent_.max = index;
  } else {
    tassert(false, "A particle can have only have up to 2 parents"
                   "(tried to add additional parent to particle '" +
                       std::string(pdg->Name()) + "')");
  }
}
//
// momentum transformations
//
// rotate from a coordate system where v moved along the z-axis, to the
// coordate system of v
// (algorithm taken from TVector3::RotateUz)
inline void particle::rotate_uz(const particle::XYZVector& v) {
  const auto uv = v.Unit();
  const double u1 = uv.X();
  const double u2 = uv.Y();
  const double u3 = uv.Z();
  double up = u1 * u1 + u2 * u2;
  if (up) {
    up = sqrt(up);
    const double px = p_.X();
    const double py = p_.Y();
    const double pz = p_.Z();
    p_ = {(u1 * u3 * px - u2 * py + u1 * up * pz) / up,
          (u2 * u3 * px + u1 * py + u2 * up * pz) / up,
          (u3 * u3 * px - px + u3 * up * pz) / up, p_.T()};
  } else if (u3 < 0) { // phi = 0, theta = pi
    p_ = {-p.X(), p.Y(), -p.Z(), p.T()};
  }
}
inline void particle::rotate_uz(const particle::XYZTVector& v) {
  rotate_uz(v.Vect());
}
inline void particle::rotate_uz(const particle& pv) {
  rotate_uz(pv.p().Vect());
}
inline void particle::boost(const particle::BetaVector& bv) {
  lorentzboost b{bv};
  boost(b);
}

//
// private utility functions
//
inline void particle::set_mass_lifetime(std::shared_ptr<TRandom> rng) {
  if (pdg->Stable()) {
    mass = pdg->Mass();
    lifetime = 0;
    else {
      mass = rng->BreitWigner(pdg->Mass(), pdg->Width());
      lifetime = rng->Exp(pdg->Lifetime());
    }
  }
}

} // namespace pcsim

#endif
