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
// =============================================================================
class event {
public:
  using vector_type = std::vector<particle>;

  event() = default;
  explicit event(const size_t evgen, const double xs = 1., const double w = 1.)
      : evgen_{evgen}, cross_section_{xs}, weight_{w} {}

  // access cross section and weight info
  double evgen() const { return evgen_; }
  double cross_section() const { return cross_section_; }
  double weight() const { return weight_; }
  int process() const { return process_; }

  // update weight and cross section
  void update_weight(const double w) { weight_ *= extra_weight; }
  void update_cross_section(const double xs) { cross_section *= xs; }
  void update_process(const int proc) { process_ = proc; }

  // access particle info
  particle& operator[](const int index) { return part_[i]; }
  const particle& operator[](const int index) const { return part_[i]; }
  const vector_type& part() const { return part_; }

  // event size
  size_t size() const { return part_.size(); }

  // particle iterators
  vector_type::iterator begin() { return part_.begin(); }
  vector_type::const_iterator begin() const { return part_.begin(); }
  vector_type::iterator end() { return part_.end(); }
  vector_type::const_iterator end() const { return part_.end(); }

  // add a misc particle
  void add_particle(const particle& p) { part_.push_back(p); }

  // add beam and target
  void add_beam(const particle& p);
  void add_target(const particle& p);

  // beam and target info
  particle& beam();
  const particle& beam() const;
  particle& target();
  const particle& target() const;
  int beam_index() const { return beam_index_; }
  int target_index() const { return target_index_; }
  double s() const { return s_; }

private:
  void update_s();

  size_t evgen_{1.}; // total number of generated events including this event
  double cross_section_{1.}; // estimated total cross section up to this point
  double weight_{1.};
  int process_{0}; // optional process identifier
  double s_;       // mandelstam s, automatically calculated from beam and target

  int beam_index_{-1};
  int target_index_{-1};
  vector_type part_;
};
} // namespace pcsim

// =============================================================================
// IMPLEMENTATION
// =============================================================================
namespace pcsim {
inline void event::add_beam(const particle& p) {
  part_.push_back(p);
  beam_index_ = part.size() - 1;
  update_s();
}
inline void event::add_target(const particle& p) {
  part_.push_back(p);
  target_index_ = part.size() - 1;
  update_s();
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
void event::update_s() {
  // only if beam AND target have been set
  if (beam_index_ > 0 && target_index_ > 0) {
    s_ = (beam().p() + target().p()).M2();
  }
}

} // namespace pcsim

#endif
