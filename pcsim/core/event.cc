#include "event.hh"
#include <pcsim/core/logger.hh>

// =============================================================================
// IMPLEMENTATION: event_out
// =============================================================================

namespace pcsim {
event_out::event_out(std::shared_ptr<TFile> f, const std::string& name)
    : file_{f} {
  LOG_INFO("event_out", "Initializing ROOT output stream");
  tassert(file_, "invalid file pointer");
  file_->cd();
  tree_ = new TTree(name.c_str(), name.c_str());
  tassert(tree_, "Failed to inialize tree " + name);
  create_branches();
}
void event_out::push(const std::vector<event>& events) {
  for (const auto e : events) {
    push(e);
  }
}
void event_out::push(const event& e) {
  evgen_ = e.evgen();
  cross_section_ = static_cast<float>(e.cross_section());
  total_cross_section_ = static_cast<float>(e.total_cross_section());
  weight_ = static_cast<float>(e.weight());
  process_ = e.process();
  s_ = static_cast<float>(e.s());
  beam_index_ = static_cast<int16_t>(e.beam_index());
  target_index_ = static_cast<int16_t>(e.target_index());

  // add the particles
  for (const auto& part : e) {
    add(part);
  }
  for (const auto& dp : e.detected()) {
    add_detected(dp);
  }

  // fill the tree
  tree_->Fill();

  // increment the index for the next event, and clear the particle buffer
  index_ += 1;
  clear();
  // that's all
}
void event_out::clear() {
  n_part_ = 0;
  part_index_.clear();
  type_.clear();
  status_.clear();
  mass_.clear();
  charge_.clear();
  p_.clear();
  vertex_.clear();
  parent_first_.clear();
  parent_second_.clear();
  daughter_begin_.clear();
  daughter_end_.clear();

  rc_n_part_ = 0;
  rc_p_.clear();
  rc_type_.clear();
  rc_index_.clear();
  rc_status_.clear();
}
void event_out::add(const particle& part) {
  n_part_ += 1;
  part_index_.push_back(static_cast<int32_t>(part.index()));
  type_.push_back(static_cast<int32_t>(part.type()));
  status_.push_back(static_cast<int16_t>(part.status()));
  mass_.push_back(static_cast<float>(part.mass()));
  charge_.push_back(static_cast<int8_t>(part.charge()));
  p_.push_back(part.p());
  vertex_.push_back(part.vertex());
  parent_first_.push_back(part.parent_first());
  parent_second_.push_back(part.parent_second());
  daughter_begin_.push_back(part.daughter_begin());
  daughter_end_.push_back(part.daughter_end());
}
// add a detected particle to the buffer
void event_out::add_detected(const detected_particle& dp) {
  rc_n_part_ += 1;
  rc_p_.push_back(dp.p());
  rc_type_.push_back(dp.generated().type<int32_t>());
  rc_index_.push_back(static_cast<int16_t>(dp.generated().index()));
  rc_status_.push_back(static_cast<int16_t>(dp.status()));
}

void event_out::create_branches() {
  tree_->Branch("index", &index_);
  tree_->Branch("evgen", &evgen_);
  tree_->Branch("cross_section", &cross_section_);
  tree_->Branch("total_cross_section", &total_cross_section_);
  tree_->Branch("weight", &weight_);
  tree_->Branch("process", &process_);
  tree_->Branch("s", &s_);
  tree_->Branch("beam_index", &beam_index_);
  tree_->Branch("target_index", &target_index_);
  tree_->Branch("n_part", &n_part_);
  tree_->Branch("part_index", &part_index_);
  tree_->Branch("type", &type_);
  tree_->Branch("status", &status_);
  tree_->Branch("mass", &mass_);
  tree_->Branch("charge", &charge_);
  tree_->Branch("p", &p_);
  tree_->Branch("vertex", &vertex_);
  tree_->Branch("parent_first", &parent_first_);
  tree_->Branch("parent_second", &parent_second_);
  tree_->Branch("daughter_begin", &daughter_begin_);
  tree_->Branch("daughter_end", &daughter_end_);
  tree_->Branch("rc_n_part", &rc_n_part_);
  tree_->Branch("rc_p", &rc_p_);
  tree_->Branch("rc_type", &rc_type_);
  tree_->Branch("rc_index", &rc_index_);
  tree_->Branch("rc_status", &rc_status_);
}

} // namespace pcsim
