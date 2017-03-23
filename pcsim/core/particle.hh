#ifndef PCSIM_CORE_PARTICLE_LOADED
#define PCSIM_CORE_PARTICLE_LOADED

#include <TLorentzVector.h>
#include <TRandom.h>
#include <TVector3.h>
#include <memory>
#include <pcsim/core/interval.hh>
#include <pcsim/core/pdg.hh>

namespace pcsim {

// particle info
class particle {
public:
  // particle status
  enum class status_code : int {
    BEAM = 11,
    SECONDARY_BEAM = 12,
    SCAT = 13,
    UNSTABLE = 21,
    UNSTABLE_SCHC = 22,
    FINAL = 31,
    DECAYED = 32,
    OTHER = 99
  };

  // constructors
  particle() = default;
  particle(const particle& rhs) = default;
  // particle with a given id
  particle(const pdg_id id)
      : type_{id}
      , pdg_{pdg_particle(id)}
      , charge_{static_cast<int>(pdg->Charge() / 3)} {}
  // particle with a given momentum 3-vector
  particle(const pdg_id id, const TVector3& p3)
      : particle{id}
      , mass_{pdg->Mass()}
      , p_{p3, sqrt(mass_ * mass_ + p_.Mag2())} {}
  // particle in a given direction with a given energy
  particle(const pdg_id id, const TVector3 direction, const double E)
      : particle{id}
      , mass_{pdg->Mass()}
      , p_{direction.Unit() * sqrt(E * E - mass_ * mass_), E} {}
  // particle with a given 4-momentum (can be off-shell)
  particle(const pdg_id id, const TLorentzVector p) {
    particle{id}, mass_{p.M()}, p_{p} {}
  }
  // RNG constructors usefull in particular for unstable particles with a mass
  // and lifetime that needs to be generated
  particle(const pdg_id id, std::shared_ptr<TRandom> rng) : particle(id) {
    set_mass_lifetime(std::move(rng));
  }
  particle(const pdg_id id, const TVector3& p3, std::shared_ptr<TRandom> rng)
      : particle{id} {
    set_mass_lifetime(std::move(rng));
    p_.SetXYZM(p.X(), p.Y(), p.Z(), mass_);
  }

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
  double lifetime() const { return lifetime_; }
  interval<int> parent() const { return parent_; }
  interval<int> daughter() const { return daughter_; }
  int daughter_begin() const { return daughter_.min; }
  int daughter_end() const { return daughter_.max; }
  int n_daughters() const { return daughter_.max - daughter_.min; }
  int parent_first() const { return parent_.min; }
  int parent_second() const { return parent_.max; }
  int n_parents() const { return parent_.max - parent_.min; }

  // get a reference to the momentum/vertex 4-vector
  TLorentzVector& p() { return p_; };
  const TLorentzVector& p() const { return p_; };
  TLorentzVector& vertex() { return vertex_; }
  const TLorentzVector& vertex() const { return vertex_; }

  // add a daughter track
  void add_daughter(const int index) {
    if (daughter_.min > index || daughter_.min < 0) {
      daughter_.min = index;
    }
    if (daughter_.max >= index || daughter_.max < 0) {
      daughter_.max = index + 1;
    }
  }
  void add_parent(const int index) {
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
  void set_parents(const interval<int> indices) { parent_ = indices; }

private:
  pdg_id type_{pdg_id::UNKNOWN};
  TParticlePDG* pdg_{nullptr};
  status_code status_{status_code::OTHER};

  int charge_{0};
  // actual mass and lifetime for this particle, can deviate from pole values
  // for unstable particles!
  double mass_{0};
  double lifetime_{0};
  TLorentzVector p_{0, 0, 0, 0};
  TLorentzVector vertex_{0, 0, 0, 0};
  // parent indices store the first (and optional second) parent of a particle.
  // -1 if not stored
  interval<int> parent_{-1, -1};
  // daughter indices are encoded from [begin, end) where begin is the
  // first index and end one past the last index. (-1, -1) if none
  interval<int> daughter_{-1, -1};
};

} // namespace pcsim

#endif
