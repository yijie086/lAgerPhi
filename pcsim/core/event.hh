#ifndef PCSIM_CORE_EVENT_LOADED
#define PCSIM_CORE_EVENT_LOADED

#include <pcsim/core/assert.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <vector>

namespace pcsim {

// =============================================================================
// event
//
// base event record that carries a list of all particles
// derive form this class for more specialized event records
//
// TODO: remove epsilon and R info, move to a gamma_p_event sub-class
// =============================================================================
class event : public generator_data {
public:
  event() = default;
  explicit event(const size_t evgen, const double xs = 1., const double w = 1.,
                 const double R = 0., const double epsilon = 0.)
      : generator_data{xs}
      , evgen_{evgen}
      , weight_{w}
      , R_{R}
      , epsilon_{epsilon} {}

  // ===========================================================================
  // EVENT INFO
  //
  // access evgen and weight info, cross section is available through the
  // generator_data base class
  size_t evgen() const { return evgen_; }
  double weight() const { return weight_; }
  double R() const { return R_; }
  double epsilon() const { return epsilon_; }
  int process() const { return process_; }
  double s() const { return s_; }

  // update weight and cross section
  void update_evgen(const size_t evgen) { evgen_ = evgen; }
  void update_weight(const double w) { weight_ *= w; }
  void update_R(const double R) { R_ = R; }
  void update_epsilon(const double epsilon) { epsilon_ = epsilon; }
  void update_process(const int proc) { process_ = proc; }

  // ===========================================================================
  // PARICLE INFO
  //
  // access particle info
  particle& operator[](const int index) { return part_[index]; }
  const particle& operator[](const int index) const { return part_[index]; }
  std::vector<particle>& part() { return part_; }
  const std::vector<particle>& part() const { return part_; }
  particle& part(const int index) { return part_[index]; }
  const particle& part(const int index) const { return part_[index]; }

  // event size
  size_t size() const { return part_.size(); }

  // particle iterators
  std::vector<particle>::iterator begin() { return part_.begin(); }
  std::vector<particle>::const_iterator begin() const { return part_.begin(); }
  std::vector<particle>::iterator end() { return part_.end(); }
  std::vector<particle>::const_iterator end() const { return part_.end(); }

  // add a misc particle. returns the index of the particle
  int add_particle(const particle& p);

  // add a daughter particle with 1 or 2 parents
  // returns the index of the daughter
  int add_daughter(particle daughter, const int parent1,
                   const int parent2 = -1);

  // add beam and target
  // returns the beam and target index
  int add_beam(const particle& p);
  int add_target(const particle& p);

  // beam and target info
  particle& beam();
  const particle& beam() const;
  particle& target();
  const particle& target() const;
  int beam_index() const { return beam_index_; }
  int target_index() const { return target_index_; }

  // ===========================================================================
  // DETECTOR INFO
  //
  // access detector info
  std::vector<detected_particle>& detected() { return detected_; }
  const std::vector<detected_particle>& detected() const { return detected_; }
  detected_particle& detected(const int index) { return detected_[index]; }
  const detected_particle& detected(const int index) const {
    return detected_[index];
  }

  int add_detected(const detected_particle dp) {
    detected_.push_back(dp);
    return detected_.size() - 1;
  }
  
private:
  void update_s();

  size_t evgen_{1}; // total number of generated events including this event
  double weight_{1.};
  double epsilon_{0.};
  double R_{0.};
  int process_{0}; // optional process identifier
  double s_; // mandelstam s, automatically calculated from beam and target

  int beam_index_{-1};
  int target_index_{-1};
  std::vector<particle> part_;
  std::vector<detected_particle> detected_;
};
} // namespace pcsim

// =============================================================================
// IMPLEMENTATION
// =============================================================================
namespace pcsim {
inline int event::add_beam(const particle& p) {
  beam_index_ = add_particle(p);
  update_s();
  return beam_index_;
}
inline int event::add_target(const particle& p) {
  target_index_ = add_particle(p);
  update_s();
  return target_index_;
}
inline int event::add_particle(const particle& p) {
  const int index = part_.size();
  part_.push_back(p);
  part_[index].update_index(index);
  return index;
}
inline int event::add_daughter(particle daughter, const int parent1,
                               const int parent2) {
  daughter.add_parent(parent1);
  daughter.add_parent(parent2);
  int index = add_particle(daughter);
  part_[parent1].add_daughter(index);
  if (parent2 >= 0) {
    part_[parent2].add_daughter(index);
  }
  return index;
}
inline particle& event::beam() {
  tassert(beam_index_ > 0, "trying to access beam data, but no beam "
                           "data present in the event.");
  return part_[beam_index_];
}
inline const particle& event::beam() const {
  tassert(beam_index_ > 0, "trying to access beam data, but no beam "
                           "data present in the event.");
  return part_[beam_index_];
}
inline particle& event::target() {
  tassert(target_index_ > 0, "trying to access target data, but no target "
                             "data present in the event.");
  return part_[target_index_];
}
inline const particle& event::target() const {
  tassert(target_index_ > 0, "trying to access target data, but no target "
                             "data present in the event.");
  return part_[target_index_];
}
inline void event::update_s() {
  // only if beam AND target have been set
  if (beam_index_ > 0 && target_index_ > 0) {
    s_ = (beam().p() + target().p()).M2();
  }
}

} // namespace pcsim

#endif
