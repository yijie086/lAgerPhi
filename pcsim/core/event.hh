#ifndef PCSIM_CORE_EVENT_LOADED
#define PCSIM_CORE_EVENT_LOADED

#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <vector>

namespace pcsim {

// base event record that carries a list of all particles
// derive form this class for more specialized event records

class event {
public:

  event() = default;
  explicit event(const double xs, const double w = 1.)
      : cross_section_{xs}, weight_{w} {}

  // access cross section and weight info
  double cross_section() const {
    return cross_section_;
  }
  double weight() const {
    return weight_;
  }

  // update weight and cross section
  void update_weight(const double w) {
    weight_ *= extra_weight;
  }
  void update_cross_section(const double xs) {
    cross_section *= xs;
  }

  // access particle info
  particle& operator[](const int index) {
    return part_[i];
  }
  const particle& operator[](const int index) const {
    return part_[i];
  }
  const std::vector<particle>& part() const {
    return part_;
  }

  // event size
  size_t size() const {
    return part_.size();
  }

  // particle iterators
  auto begin() {
    return part_.begin();
  }
  auto end() {
    return part_.end();
  }

  // add a misc particle
  void add_particle(const particle& p) { part_.push_back(p); }

  // add beam and target
  void add_beam(const particle& p) {
    part_.push_back(p);
    beam_index_ = part.size() - 1;
  }
  void add_target(const particle& p) {
    part_.push_back(p);
    target_index_ = part.size() - 1;
  }

  // beam and target info
  particle& beam() {
    tassert(beam_index_ > 0, "trying to access beam data, but no beam "
                             "data present in the event.");
    return part_[beam_index_];
  }
  const particle& beam() const {
    tassert(beam_index_ > 0, "trying to access beam data, but no beam "
                             "data present in the event.");
    return part_[beam_index_];
  }
  int beam_index() const { return beam_index_; }
  particle& target() {
    tassert(target_index_ > 0, "trying to access target data, but no target "
                               "data present in the event.");
    return part_[target_index_];
  }
  const particle& target() const {
    tassert(target_index_ > 0, "trying to access target data, but no target "
                               "data present in the event.");
    return part_[target_index_];
  }
  int target_index() const { return target_index_; }

private:
  cross_section_{1.};
  weight_{1.};
  int beam_index_{-1};
  int target_index_{-1};
  std::vector<particle> part_;
};

}

#endif
