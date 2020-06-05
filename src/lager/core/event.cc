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

#include "event.hh"
#include <cstring>
#include <lager/core/logger.hh>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

// =============================================================================
// IMPLEMENTATION: event_out
// =============================================================================

namespace lager {
event_out::event_out(std::shared_ptr<TFile> f,
                     std::unique_ptr<HepMC3::WriterAscii> ohepmc,
                     std::unique_ptr<std::ofstream> ogemc,
                     std::unique_ptr<std::ofstream> osimc,
                     const std::string& name)
    : file_{f}
    , ohepmc_{std::move(ohepmc)}
    , ogemc_{std::move(ogemc)}
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
  ibeam_index_ = static_cast<int16_t>(e.ibeam_index());
  tbeam_index_ = static_cast<int16_t>(e.tbeam_index());

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

  // write HEPMC record if wanted
  if (ohepmc_) {
    write_hepmc(e);
  }

  // write GEMC record if wanted
  if (ogemc_) {
    write_gemc(e);
  }
  // write SIMC record if wanted
  if (osimc_) {
    write_simc(e);
  }

  // that's all
}

void event_out::write_hepmc(const event& e) {
  // event index
  static size_t index = 0;
  // create event
  HepMC3::GenEvent hevt(HepMC3::Units::GEV, HepMC3::Units::CM);
  // get vector of HepMC particles corresponding to our particles
  std::vector<HepMC3::GenParticlePtr> hepmc_part;
  std::transform(e.part().begin(), e.part().end(),
                 std::back_inserter(hepmc_part), [](const particle& part) {
                   // final state: 1
                   // decayed: 2
                   // documentation: 3
                   const int status =
                       part.final_state()
                           ? 1
                           : (part.decayed() ? 2
                                             : (!part.documentation() ? 3 : 0));
                   return HepMC3::GenParticlePtr(new HepMC3::GenParticle(
                       HepMC3::FourVector(part.p().X(), part.p().Y(),
                                          part.p().Z(), part.p().E()),
                       part.type<int>(), status));
                 });
  // record of "processed" particles. A particle is processed once it is added
  // as incoming leg of a vertex
  std::vector<size_t> finished;
  // loop over all particles, find and create vertices and build the HepMC event
  for (size_t i = 0; i < e.part().size(); ++i) {
    const auto& part = e.part(i);
    // 0. Only create vertices from "incoming" lines
    if (e.part(i).n_daughters() == 0) {
      continue;
    }
    const auto& first_daughter = e.part(part.daughter_begin());
    // 1. Check if event is already "processed"
    if (std::any_of(finished.begin(), finished.end(),
                    [i](int j) { return (j == i); })) {
      continue;
    }
    // 2. OK, let's create a new vertex for this particle
    //    In principle the vertex member of a particle is the start vertex,
    //    so we get the relevant vertex from the first daughter particle instead
    auto raw_vertex = first_daughter.vertex();
    auto vx = HepMC3::GenVertexPtr(new HepMC3::GenVertex(HepMC3::FourVector(
        raw_vertex.X(), raw_vertex.Y(), raw_vertex.Z(), raw_vertex.T())));
    // 3. Attach incoming lines to this vertex and mark them as "finished"
    //    Use the first daughter to get the full list of incoming lines
    for (int iin = first_daughter.parent_first();
         iin < first_daughter.parent_second(); ++iin) {
      vx->add_particle_in(hepmc_part[iin]);
      finished.push_back(iin);
    }
    // 4. attach outgoing line to this vertex
    for (int iout = part.daughter_begin(); iout < part.daughter_end(); ++iout) {
      vx->add_particle_out(hepmc_part[iout]);
    }
    // 5. Add vertex to event
    hevt.add_vertex(vx);
  }
  // Now we are ready to write out the event
  ohepmc_->write_event(hevt);
}

void event_out::write_gemc(const event& e) {
  char buf[2048];
  // write the first line of the event record
  snprintf(buf, 2048, "%5zu %5i %5i %5f %5f %5i %5f %5i %5i %8.6e\n",
           e.count_final_state(), 1, 1, 0., 0., 11, e.ibeam().energy(), 0,
           e.process(), e.total_cross_section());
  *ogemc_ << buf;
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
    *ogemc_ << buf;
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
             hms_track->p().X(), hms_track->p().Y(), hms_track->p().Z(),
             hms_track->energy(), hms_track->vertex().Z(), shms_track->p().X(),
             shms_track->p().Y(), shms_track->p().Z(), shms_track->energy(),
             shms_track->vertex().Z(), 1.);
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
  tree_->Branch("ibeam_index", &ibeam_index_);
  tree_->Branch("tbeam_index", &tbeam_index_);
  tree_->Branch("n_part", &n_part_);
  tree_->Branch("particles", &parts_);
  tree_->Branch("rc_n_part", &rc_n_part_);
  tree_->Branch("rc_particles", &rc_parts_);
}

} // namespace lager
