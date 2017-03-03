#ifndef PCSIM_CORE_PARTICLE_LOADED
#define PCSIM_CORE_PARTICLE_LOADED

#include <TLorentzVector.h>
#include <TRandom.h>
#include <TVector3.h>
#include <memory>
#include <pcsim/core/interval.hh>
#include <pcsim/core/pdg.hh>

namespace pcsim {

// basic particle info 
struct particle_info {
  TParticlePDG* pdg; // pdg particle info
  int32_t info;      // LUND particle code
  int32_t charge;    // particle charge

  particle_info() = default;
  explicit particle_info(const pdg_id id)
      : pdg{pdg_particle(id)}
      , info{static_cast<int32_t>(id)}
      , charge{static_cast<int32_t>(pdg->Charge() / 3)} {}
  particle_info(const particle_info& rhs) = default;

  std::string name() const { return pdg->GetName(); }
};

// particle, will initialize mass and lifetime for unstable particles if
// a random generator is provided in the constructor
struct particle : particle_info {
  double mass = 0;                      // mass in GeV/c^2
  double lifetime = 0;                  // lifetime in cm/c
  TLorentzVector mom = {0, 0, 0, 0};    // momentum in GeV/c
  TLorentzVector vertex = {0, 0, 0, 0}; // vertex in cm

  // normal constructors
  particle() = default;
  explicit particle(const pdg_id id)
      : particle_info{id}, mass{pdg->Mass()}, mom{0, 0, 0, mass} {}
  particle(const pdg_id id, const TVector3& p)
      : particle_info{id}
      , mass{pdg->Mass()}
      , mom{p, sqrt(mass * mass + p.Mag2())} {}
  particle(const pdg_id id, const TVector3& dir, const double E)
      : particle_info{id}
      , mass{pdg->Mass()}
      , mom{dir.Unit() * sqrt(E * E - mass * mass), E} {}
  particle(const particle& rhs)
      : particle_info{static_cast<const particle_info&>(rhs)}
      , mass{rhs.mass}
      , lifetime{rhs.lifetime}
      , mom{rhs.mom}
      , vertex{rhs.vertex} {}
  explicit particle(const particle_info& rhs)
      : particle_info{rhs}, mass{pdg->Mass()} {}
  // off-shell particle
  particle(const pdg_id id, const TLorentzVector& pp)
      : particle_info{id}, mass{pp.M()}, mom{pp} {}
  // RNG constructors that generate mass and/or lifetime
  particle(const pdg_id id, std::shared_ptr<TRandom> rng) : particle_info{id} {
    set_mass_lifetime(std::move(rng));
  }
  particle(const pdg_id id, const TVector3& p, std::shared_ptr<TRandom> rng)
      : particle_info{id} {
    set_mass_lifetime(std::move(rng));
    mom.SetXYZM(p.X(), p.Y(), p.Z(), mass);
  }
  //particle(const pdg_id id, const TVector3& dir, const double E,
  //         std::shared_ptr<TRandom> rng)
  //    : particle_info{id} {
  //  set_mass_lifetime(std::move(rng));
  //  mom = {dir.Unit() * sqrt(E * E - mass * mass), E};
  //}
  //particle(const particle_info& rhs, std::shared_ptr<TRandom> rng)
      //: particle_info{rhs} {
    //set_mass_lifetime(std::move(rng));
  //}

  void set_lifetime(std::shared_ptr<TRandom> rng) {
    if (pdg->Width() > 0) {
      lifetime = rng->Exp(pdg->Lifetime()) * TMath::Ccgs();
    } else {
      lifetime = 0;
    }
  }
  void set_mass_lifetime(std::shared_ptr<TRandom> rng) {
    if (pdg->Width() > 0) {
      mass = rng->BreitWigner(pdg->Mass(), pdg->Width());
      lifetime = rng->Exp(pdg->Lifetime()) * TMath::Ccgs();
    } else {
      mass = pdg->Mass();
      lifetime = 0;
    }
  }
};

// full MC particle info, based on the LUND event record
struct mc_particle : particle {
  enum status_code { INITIAL, VIRTUAL, UNSTABLE, FINAL, SCAT, OTHER };
  status_code status = OTHER; // particle status
  // parent/daughter indeces are encoded from [begin, end] where begin is the
  // first index and end one past the last index. (-1, -1) if none
  interval<int> parent = {-1, -1};
  interval<int> daughter = {-1, -1};

  // normal constructors
  explicit mc_particle(const pdg_id id) : particle{id} {}
  mc_particle(const pdg_id id, const status_code status)
      : particle{id}, status{status} {}
  mc_particle(const mc_particle& rhs)
      : particle{static_cast<const particle&>(rhs)}
      , status{rhs.status}
      , parent{rhs.parent}
      , daughter{rhs.daughter} {}
  explicit mc_particle(const particle& rhs) : particle{rhs} {}
  mc_particle(const particle& rhs, const status_code status)
      : particle{rhs}, status{status} {}

  void add_daughter(int index) {
    if (daughter.min > index || daughter.min < 0) {
      daughter.min = index;
    }
    if (daughter.max >= index || daughter.max < 0) {
      daughter.max = index + 1;
    }
  }
  void set_parents(interval<int> p) {
    parent = p;
  }
};

} // namespace pcsim

#endif
