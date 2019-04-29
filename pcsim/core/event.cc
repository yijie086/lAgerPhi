#include "event.hh"
#include <cstring>
#include <pcsim/core/logger.hh>

// =============================================================================
// IMPLEMENTATION: event_out
// =============================================================================

namespace pcsim {
event_out::event_out(std::shared_ptr<TFile> f,
                     std::unique_ptr<std::ofstream> olund,
                     std::unique_ptr<std::ofstream> osimc,
                     const std::string& name)
    : file_{f}
    , olund_{std::move(olund)}
    , osimc_{std::move(osimc)}
    , parts_{"TParticle", PARTICLE_BUFFER_SIZE}
    , rc_parts_{"TParticle", PARTICLE_BUFFER_SIZE} {
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
  process_cross_section_ = static_cast<float>(e.process_cross_section());
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

  // write LUND record if wanted
  if (olund_) {
    write_lund(e);
  }
  // write SIMC record if wanted
  if (osimc_) {
    write_simc(e);
  }

  // that's all
}

void event_out::write_lund(const event& e) {
  char buf[2048];
  // write the first line of the event record
  snprintf(buf, 2048, "%5zu %5i %5i %5f %5f %5i %5f %5i %5i %8.6e\n",
           e.count_final_state(), 1, 1, 0., 0., 11, e.beam().energy(), 0,
           e.process(), e.total_cross_section());
  *olund_ << buf;
  for (const auto& part : e) {
    if (!part.final_state()) {
      continue;
    }
    snprintf(buf, 2048,
             "%5i %5f  %5i %5i %5i %5i %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e "
             "%8.6e %8.6e\n",
             part.index(), part.lifetime(), 1, part.type<int>(),
             part.parent_first(), part.daughter_begin(), part.p().x(),
             part.p().y(), part.p().z(), part.energy(), part.mass(),
             part.vertex().x(), part.vertex().y(), part.vertex().z());
    *olund_ << buf;
  }
}
void event_out::write_simc(const event& e) {
  char buf[1024];
  // assume HMS is detector ID 1 and SHMS is detector ID 2
  auto hms_track =
      std::find_if(e.detected().begin(), e.detected().end(),
                   [](const auto& part) { return part.status() == 1; });
  auto shms_track =
      std::find_if(e.detected().begin(), e.detected().end(),
                   [](const auto& part) { return part.status() == 2; });
  auto end = e.detected().end();
  if (hms_track != end && shms_track != end) {
    snprintf(buf, 1024,
             "%16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e "
             "%16.10e %16.10e %16.10e\n",
             hms_track->p().X(), hms_track->p().Y(),
             hms_track->p().Z(), hms_track->energy(),
             hms_track->vertex().Z(), shms_track->p().X(),
             shms_track->p().Y(), shms_track->p().Z(),
             shms_track->energy(), shms_track->vertex().Z(), 1.);
    *osimc_ << buf;
  }
}

void event_out::clear() {
  n_part_ = 0;
  rc_n_part_ = 0;
  parts_.Clear();
  rc_parts_.Clear();
}
void event_out::add(const particle& part) {
  auto pbuf = new (parts_[n_part_]) TParticle(
      static_cast<int32_t>(part.type()), static_cast<int32_t>(part.status()),
      part.parent_first(), part.parent_second(), part.daughter_begin(),
      part.daughter_end(), part.p().X(), part.p().Y(), part.p().Z(),
      part.p().E(), part.vertex().X(), part.vertex().Y(), part.vertex().Z(),
      part.vertex().T());
  // set the weight of non final state particles to zero
  pbuf->SetWeight(part.final_state() ? 1. : 0.);
  // increment our particle counter
  n_part_ += 1;
}
// add a detected particle to the buffer
void event_out::add_detected(const detected_particle& dp) {
  auto pbuf = new (rc_parts_[rc_n_part_])
      TParticle(static_cast<int32_t>(dp.generated().type<int32_t>()),
                dp.status(), dp.generated().index(), 0, 0, 0, dp.p().X(),
                dp.p().Y(), dp.p().Z(), dp.p().E(), dp.vertex().X(),
                dp.vertex().Y(), dp.vertex().Z(), dp.vertex().T());
  // increment our detected particle counter
  rc_n_part_ += 1;
}

void event_out::create_branches() {
  tree_->Branch("index", &index_);
  tree_->Branch("evgen", &evgen_);
  tree_->Branch("cross_section", &cross_section_);
  tree_->Branch("total_cross_section", &total_cross_section_);
  tree_->Branch("process_cross_section", &process_cross_section_);
  tree_->Branch("weight", &weight_);
  tree_->Branch("process", &process_);
  tree_->Branch("s", &s_);
  tree_->Branch("beam_index", &beam_index_);
  tree_->Branch("target_index", &target_index_);
  tree_->Branch("n_part", &n_part_);
  tree_->Branch("particles", &parts_);
  tree_->Branch("rc_n_part", &rc_n_part_);
  tree_->Branch("rc_particles", &rc_parts_);
}

} // namespace pcsim
