// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2020 Sylvester Joosten <sjoosten@anl.gov>
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

#ifndef LAGER_CORE_PARTICLE_LOADED
#define LAGER_CORE_PARTICLE_LOADED

#include <Math/LorentzRotation.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TRandom.h>
#include <lager/core/assert.hh>
#include <lager/core/interval.hh>
#include <lager/core/pdg.hh>
#include <memory>

namespace lager {

// =============================================================================
// particle
//
// cornerstone class containing particle info
// =============================================================================
class particle {
public:
  using BetaVector = ROOT::Math::XYZTVector::BetaVector;
  using XYZTVector = ROOT::Math::XYZTVector;
  using Boost = ROOT::Math::Boost;
  using XYZVector = ROOT::Math::XYZVector;
  using Polar3DVector = ROOT::Math::Polar3DVector;

  // particle status
  enum class status_code : int {
    BEAM = 11,
    SECONDARY_BEAM = 12,
    SCAT = 13,
    RECOIL = 14,
    SPECTATOR = 15,
    DECAYED = 21,
    DECAYED_SCHC = 22,
    DECAYED_RADCOR_ONLY = 23,
    FINAL = 30,
    UNSTABLE = 31,
    UNSTABLE_SCHC = 32,
    UNSTABLE_RADCOR_ONLY = 33,
    INFO = 40,
    INFO_PARENT_CM = 41,
    OTHER = 99
  };

  // constructors
  //
  // default
  particle(const particle&) = default;
  particle& operator=(const particle&) = default;
  // particle with a given id or name or int
  // All other constructors are templates that will resolve to one
  // of these types
  particle(const pdg_id id = pdg_id::unknown,
           const status_code status = status_code::FINAL);
  particle(const std::string& name,
           const status_code status = status_code::FINAL);
  particle(const int32_t id, const status_code status = status_code::FINAL);
  // particle with a given momentum 3-vector
  template <class PdgId>
  particle(const PdgId& id, const XYZVector& p3,
           const status_code status = status_code::FINAL);

  // particle in a given direction with a given energy
  template <class PdgId>
  particle(const PdgId& id, const XYZVector& direction, const double E,
           const status_code status = status_code::FINAL);

  // particle with a given 4-momentum (can be off-shell)
  template <class PdgId>
  particle(const PdgId& id, const XYZTVector& p,
           const status_code status = status_code::FINAL);

  // RNG constructors usefull in particular for unstable particles with a mass
  // and lifetime that needs to be generated
  // status will be auto-set to UNSTABLE for unstable particles, and FINAL for
  // stable particles
  template <class PdgId>
  particle(const PdgId& id, std::shared_ptr<TRandom> rng);
  template <class PdgId>
  particle(const PdgId& id, const XYZVector& p3, std::shared_ptr<TRandom> rng);

  // particle ID code
  template <class Integer = pdg_id> Integer type() const {
    return static_cast<Integer>(type_);
  }

  // particle index, -1 if not a part of a structured event
  int index() const { return index_; }
  void update_index(const int i) { index_ = i; };

  // full PDG info from the TDatabasePDG
  TParticlePDG* pdg() const { return pdg_; }

  // particle status
  template <class Integer = status_code> Integer status() const {
    return static_cast<Integer>(status_);
  }
  void update_status(const status_code ns) { status_ = ns; }
  bool stable() const {
    return (status_ != status_code::UNSTABLE &&
            status_ != status_code::UNSTABLE_SCHC &&
            status_ != status_code::UNSTABLE_RADCOR_ONLY &&
            status_ != status_code::DECAYED &&
            status_ != status_code::DECAYED_SCHC &&
            status_ != status_code::DECAYED_RADCOR_ONLY);
  }
  bool decayed() const {
    return (status_ == status_code::DECAYED ||
            status_ == status_code::DECAYED_SCHC ||
            status_ == status_code::DECAYED_RADCOR_ONLY);
  }

  // final state particles (labeled FINAL or SCAT, RECOIL or SPECTATOR)
  // Also includes undecayed unstable particles
  bool final_state() const {
    return (status_ == status_code::FINAL || status_ == status_code::SCAT ||
            status_ == status_code::SPECTATOR ||
            status_ == status_code::RECOIL ||
            status_ == status_code::UNSTABLE ||
            status_ == status_code::UNSTABLE_SCHC ||
            status_ == status_code::UNSTABLE_RADCOR_ONLY);
  }

  bool documentation() const {
    return (status_ == status_code::INFO) ||
           (status_ == status_code::INFO_PARENT_CM);
  }

  // particle properties
  //
  // charge
  int charge() const { return charge_; }
  // mass and mass squared, will be taken from 4-vector of RNG if requested
  double mass() const { return mass_; }
  double mass2() const { return mass_ * mass_; }
  // pole mass (differs from mass for unstable particles
  double pole_mass() const { return pdg_->Mass(); }
  // width for unstable particles
  double width() const { return width_; }
  // actual generated lifetime for unstable particles
  double lifetime() const { return lifetime_; }
  // parent or parents
  interval<int> parent() const { return parent_; }
  // first and last+1 index of the daughters
  interval<int> daughter() const { return daughter_; }
  // other parent-daughter accessors
  int daughter_begin() const { return daughter_.min; }
  int daughter_end() const { return daughter_.max; }
  int n_daughters() const { return daughter_.max - daughter_.min; }
  int parent_first() const { return parent_.min; }
  int parent_second() const { return parent_.max; }
  int n_parents() const;
  // momentum and energy
  double momentum() const {
    return sqrt(p_.X() * p_.X() + p_.Y() * p_.Y() + p_.Z() * p_.Z());
  }
  double energy() const { return p_.E(); }
  double theta() const { return p_.theta(); }
  double phi() const { return p_.phi(); }
  // name
  std::string name() const { return pdg_->GetName(); }

  // get a reference to the momentum/vertex 4-vector
  XYZTVector& p() { return p_; };
  const XYZTVector& p() const { return p_; };
  XYZTVector& vertex() { return vertex_; }
  const XYZTVector& vertex() const { return vertex_; }

  // transformations etc.
  void boost(const BetaVector& bv);
  void boost(const Boost& b) { p_ = b * p_; }
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
    return {v.X(), v.Y(), v.Z(), t};
  }

  int index_{-1};
  pdg_id type_{pdg_id::unknown};
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
// =============================================================================
// DETECTED PARTICLE
//
// small utility class for detected particles
// =============================================================================
class detected_particle {
public:
  using XYZTVector = particle::XYZTVector;

  // constructors
  //
  // default
  detected_particle(const detected_particle&) = default;
  detected_particle& operator=(const detected_particle&) = default;
  // from a particle with a given index
  detected_particle(const particle& part, const XYZTVector& momentum,
                    const XYZTVector vertex, const int status = 0)
      : status_{status}, generated_{&part}, p_{momentum}, vertex_{vertex} {}
  detected_particle(const particle& part, const XYZTVector& momentum,
                    const int status = 0)
      : status_{status}
      , generated_{&part}
      , p_{momentum}
      , vertex_{part.vertex()} {}
  detected_particle(const particle& part, const int status = 0)
      : status_{status}
      , generated_{&part}
      , p_{part.p()}
      , vertex_{part.vertex()} {}

  // particle status
  int status() const { return status_; }
  void update_status(const int ns) { status_ = ns; }

  const particle& generated() const {
    tassert(generated_, "Associated generated particle to this detected "
                        "particle is a nullptr");
    return *generated_;
  }

  // detected 4-momentum
  const XYZTVector& p() const { return p_; }
  // detected vertex
  const XYZTVector& vertex() const { return vertex_; }

  // mass
  double mass() const { return p_.M(); }
  double mass2() const { return p_.M2(); }
  // momentum and energy
  double momentum() const { return sqrt(p_.Vect().Mag2()); }
  double energy() const { return p_.E(); }

private:
  int status_{0};
  const particle* generated_{nullptr};
  XYZTVector p_;
  XYZTVector vertex_;
};

} // namespace lager

// =============================================================================
// Particle implementation
//
// This is a work-horse class called, and therefore all its members are defined
// inline
// =============================================================================
namespace lager {

//
// CONSTRUCTORS
//
// particle with a given id
inline particle::particle(const pdg_id id, const status_code status)
    : type_{id}
    , pdg_{pdg_particle(id)}
    , charge_{static_cast<int>(pdg_->Charge() / 3)}
    , width_{pdg_->Width()}
    , status_{status}
    , mass_{pdg_->Mass()}
    , p_{0, 0, 0, mass_} {}
// particle with a given name
inline particle::particle(const std::string& name, const status_code status)
    : pdg_{pdg_particle(name)}
    , charge_{static_cast<int>(pdg_->Charge() / 3)}
    , width_{pdg_->Width()}
    , status_{status}
    , mass_{pdg_->Mass()}
    , p_{0, 0, 0, mass_} {
  type_ = static_cast<pdg_id>(pdg_->PdgCode());
}
inline particle::particle(const int32_t id, const status_code status)
    : particle(static_cast<pdg_id>(id), status) {}
// particle with a given momentum 3-vector
template <class PdgId>
particle::particle(const PdgId& id, const XYZVector& p3,
                   const status_code status)
    : particle{id, status} {
  p_ = make_4vector(p3, sqrt(mass_ * mass_ + p3.Mag2()));
}
// particle in a given direction with a given energy
template <class PdgId>
particle::particle(const PdgId& id, const XYZVector& direction, const double E,
                   const status_code status)
    : particle{id, status} {
  p_ = make_4vector(direction.Unit() * sqrt(E * E - mass_ * mass_), E);
}
// particle with a given 4-momentum (can be off-shell)
template <class PdgId>
particle::particle(const PdgId& id, const XYZTVector& p,
                   const status_code status)
    : particle{id, status} {
  mass_ = p.M();
  p_ = p;
}
// RNG constructors usefull in particular for unstable particles with a mass
// and lifetime that needs to be generated
template <class PdgId>
particle::particle(const PdgId& id, std::shared_ptr<TRandom> rng)
    : particle(id) {
  set_mass_lifetime(std::move(rng));
}
template <class PdgId>
particle::particle(const PdgId& id, const XYZVector& p3,
                   std::shared_ptr<TRandom> rng)
    : particle{id} {
  set_mass_lifetime(std::move(rng));
  p_.SetXYZT(p3.X(), p3.Y(), p3.Z(), sqrt(mass_ * mass_ + p3.Mag2()));
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
                       name() + "')");
  }
}

inline int particle::n_parents() const {
  if (parent_.min < 0) {
    return 0;
  } else if (parent_.max < 0) {
    return 1;
  } else {
    return 2;
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
    p_ = {-p_.X(), p_.Y(), -p_.Z(), p_.T()};
  }
}
inline void particle::rotate_uz(const particle::XYZTVector& v) {
  rotate_uz(v.Vect());
}
inline void particle::rotate_uz(const particle& pv) {
  rotate_uz(pv.p().Vect());
}
inline void particle::boost(const particle::BetaVector& bv) {
  Boost b{bv};
  boost(b);
}

//
// private utility functions
//
inline void particle::set_mass_lifetime(std::shared_ptr<TRandom> rng) {
  if (pdg_->Stable()) {
    mass_ = pdg_->Mass();
    lifetime_ = 0;
    status_ = status_code::FINAL;
  } else {
    mass_ = rng->BreitWigner(pdg_->Mass(), pdg_->Width());
    lifetime_ = rng->Exp(pdg_->Lifetime());
    status_ = status_code::UNSTABLE;
  }
}

} // namespace lager

#endif
