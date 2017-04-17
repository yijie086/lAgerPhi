#include "event_out.hh"
#include <pcsim/core/assert.hh>
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
void event_out::push(const event& e) {
  evgen_ = e.evgen();
  cross_section_ = static_cast<float>(e.cross_section());
  epsilon_ = static_cast<float>(e.epsilon());
  R_ = static_cast<float>(e.R());
  weight_ = static_cast<float>(e.weight());
  process_ = e.process();
  s_ = static_cast<float>(e.s());
  beam_index_ = static_cast<int16_t>(e.beam_index());
  target_index_ = static_cast<int16_t>(e.target_index());

  // add the particles
  for (const auto& part : e) {
    add(part);
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
}
void event_out::add(const particle& part) {
  n_part_ += 1;
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
void event_out::create_branches() {
  tree_->Branch("index", &index_);
  tree_->Branch("evgen", &evgen_);
  tree_->Branch("cross_section", &cross_section_);
  tree_->Branch("epsilon", &epsilon_);
  tree_->Branch("R", &R_);
  tree_->Branch("weight", &weight_);
  tree_->Branch("process", &process_);
  tree_->Branch("s", &s_);
  tree_->Branch("beam_index", &beam_index_);
  tree_->Branch("target_index", &target_index_);
  tree_->Branch("n_part", &n_part_);
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
}

} // namespace pcsim
