// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
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

#ifndef LAGER_CORE_EVENT_LOADED
#define LAGER_CORE_EVENT_LOADED

#include <lager/core/assert.hh>
#include <lager/core/event.hh>
#include <lager/core/generator.hh>
#include <lager/core/particle.hh>

#include <TClonesArray.h>
#include <TFile.h>
#include <TParticle.h>
#include <TTree.h>

#include <HepMC3/WriterAscii.h>

#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace lager {

// =============================================================================
// event
//
// base event record that carries a list of all particles
// derive form this class for more specialized event records
// =============================================================================
class event : public generator_data {
public:
  event(const event&) = default;
  event& operator=(const event&) = default;
  explicit event(const double xs = 1., const double w = 1.)
      : generator_data{xs}, weight_{w} {}

  // ===========================================================================
  // EVENT INFO
  //
  // access evgen and weight info, cross section is available through the
  // generator_data base class
  size_t evgen() const { return evgen_; }
  double total_cross_section() const { return total_cross_section_; }
  double process_cross_section() const { return process_cross_section_; }
  double weight() const { return weight_; }
  int process() const { return process_; }

  // event CM energy (from beam/target)
  double s() const;

  // update MC info (evgen, total cross section and process cross
  // section, also apply jacobians if necessary to go to a more sane
  // differential form, in case we generated in some clever unintuitive phase
  // space)
  void update_stat(const size_t evgen, const double txs);
  void update_process(const int proc);
  // update weight and process number
  void update_weight(const double w) { weight_ *= w; }
  void reset_weight(const double w = 1.) { weight_ = w; }

  // count the number of tracks with certain properties
  size_t count_status(const particle::status_code stat) const {
    return std::count_if(part_.begin(), part_.end(),
                         [=](const particle& p) { return p.status() == stat; });
  }
  size_t count_stable() const {
    return std::count_if(part_.begin(), part_.end(),
                         [](const particle& p) { return p.stable(); });
  }
  size_t count_final_state() const {
    return std::count_if(part_.begin(), part_.end(),
                         [](const particle& p) { return p.final_state(); });
  }

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
  int add_particle(const std::pair<particle, particle>& p);

  // add a daughter particle with 1 or 2 parents
  // returns the index of the daughter
  int add_daughter(particle daughter, const int parent1,
                   const int parent2 = -1);
  int add_daughter(const std::pair<particle, particle>& daughters,
                   const int parent1, const int parent2 = -1);

  // add incoming and target beams
  // returns the beam and target index
  int add_ibeam(const particle& p);
  int add_tbeam(const particle& p);

  // beam and target info
  particle& ibeam();
  const particle& ibeam() const;
  particle& tbeam();
  const particle& tbeam() const;
  int ibeam_index() const { return ibeam_index_; }
  int tbeam_index() const { return tbeam_index_; }

  // ===========================================================================
  // DETECTOR INFO
  //
  // access detector info
  std::vector<detected_particle>& detected() { return detected_; }
  const std::vector<detected_particle>& detected() const { return detected_; }
  detected_particle& detected(const int index) { return detected_[index]; }
  const detected_particle& detected(const int index) const;

  int add_detected(const detected_particle& dp);

private:
  size_t evgen_{1}; // total number of generated events including this event
  double total_cross_section_{0.};   // estimated total integrated cross section
  double process_cross_section_{0.}; // differential process cross section
  double weight_{1.};
  int process_{0}; // optional process identifier

  double s_; // mandelstam s, automatically calculated from beam and target

  int ibeam_index_{-1};
  int tbeam_index_{-1};

  //
  std::vector<particle> part_;
  std::vector<detected_particle> detected_;
};
} // namespace lager

// =============================================================================
// event_out
//
// Output for the base event record.
//
// Derive from this class for more specialized event records.
// (where the specialized event record derives from the main event class)
//    * Make sure to define your own push(your_event_type) method, and call the
//      parent push(parent_event_type) from within the method.
//    * you are responsible to create the necessary branches for your custom
//      event type, the main event branches are added by this base class
// =============================================================================

// TODO needs refactoring of the output plugins
//      --> migrate to config-based approach rather
//          than hardcoded compontents
namespace lager {
class event_out {
public:
  constexpr static const int32_t PARTICLE_BUFFER_SIZE{1000};

  event_out(std::shared_ptr<TFile> f,
            std::unique_ptr<HepMC3::WriterAscii> ohepmc,
            std::unique_ptr<std::ofstream> ogemc,
            std::unique_ptr<std::ofstream> osimc, const std::string& name);
  ~event_out() { tree_->AutoSave(); }

  // no implicit default constructors
  event_out() = delete;
  event_out(const event_out&) = delete;
  event_out& operator=(const event_out&) = delete;

  // add a event(s) to the event buffer, and flush the buffer to the tree
  void push(const event& e);
  void push(const std::vector<event>& e);

  TTree* tree() { return tree_; }

private:
  void write_hepmc(const event& e);
  void write_gemc(const event& e);
  void write_simc(const event& e);

  // clear particle portion of the event buffer
  void clear();
  // add a particle to the event buffer
  void add(const particle& part);
  // add a detected particle to the buffer
  void add_detected(const detected_particle& dp);
  // create branches
  void create_branches();

private:
  // file and tree
  std::shared_ptr<TFile> file_;
  TTree* tree_; // raw pointer because the TFile will have ownership of the tree
  std::unique_ptr<HepMC3::WriterAscii> ohepmc_; // HEPMC output stream
  std::unique_ptr<std::ofstream> ogemc_;       // GEMC output stream
  std::unique_ptr<std::ofstream> osimc_;       // SIMC output stream

  // event data
  int32_t index_{0};
  int32_t evgen_;
  float cross_section_;
  float total_cross_section_;
  float process_cross_section_;
  float weight_;
  int32_t process_;
  float s_;
  int16_t ibeam_index_;
  int16_t tbeam_index_;

  // particle data
  int16_t n_part_{0};
  TClonesArray parts_;
  int16_t rc_n_part_{0};
  TClonesArray rc_parts_;
};
} // namespace lager

// =============================================================================
// EVENT IMPLEMENTATION
// =============================================================================
namespace lager {
inline int event::add_particle(const particle& p) {
  const int index = part_.size();
  part_.push_back(p);
  part_[index].update_index(index);
  return index;
}
inline int event::add_particle(const std::pair<particle, particle>& p) {
  int first = add_particle(p.first);
  add_particle(p.second);
  return first;
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
inline int event::add_daughter(const std::pair<particle, particle>& daughters,
                               const int parent1, const int parent2) {
  int first = add_daughter(daughters.first, parent1, parent2);
  add_daughter(daughters.second, parent1, parent2);
  return first;
}

inline void event::update_stat(const size_t evgen, const double txs) {
  evgen_ = evgen;
  total_cross_section_ = txs;
  update_cross_section(jacobian());
  process_cross_section_ = 0.; // not used anymore
}
inline void event::update_process(const int proc) { process_ = proc; }
inline const detected_particle& event::detected(const int index) const {
  return detected_[index];
}

inline int event::add_detected(const detected_particle& dp) {
  detected_.push_back(dp);
  return detected_.size() - 1;
}
inline int event::add_ibeam(const particle& p) {
  ibeam_index_ = add_particle(p);
  return ibeam_index_;
}
inline int event::add_tbeam(const particle& p) {
  tbeam_index_ = add_particle(p);
  return tbeam_index_;
}
inline particle& event::ibeam() {
  tassert(ibeam_index_ >= 0, "trying to access beam data, but no beam "
                             "data present in the event.");
  return part_[ibeam_index_];
}
inline const particle& event::ibeam() const {
  tassert(ibeam_index_ >= 0, "trying to access beam data, but no beam "
                             "data present in the event.");
  return part_[ibeam_index_];
}
inline particle& event::tbeam() {
  tassert(tbeam_index_ >= 0, "trying to access target data, but no target "
                             "data present in the event.");
  return part_[tbeam_index_];
}
inline const particle& event::tbeam() const {
  tassert(tbeam_index_ >= 0, "trying to access target data, but no target "
                             "data present in the event.");
  return part_[tbeam_index_];
}
inline double event::s() const {
  // only if beam AND target have been set
  if (ibeam_index_ < 0 || tbeam_index_ < 0) {
    return 0.;
  }
  return (ibeam().p() + tbeam().p()).M2();
}

} // namespace lager

#endif
