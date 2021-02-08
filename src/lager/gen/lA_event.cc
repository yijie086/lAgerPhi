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

#include "lA_event.hh"

namespace lager {

// =============================================================================
// IMPLEMENTATION: event_out
// =============================================================================
lA_out::lA_out(std::shared_ptr<TFile> f,
               std::unique_ptr<HepMC3::WriterAscii> ohepmc,
               std::unique_ptr<std::ofstream> ogemc,
               std::unique_ptr<std::ofstream> osimc, const std::string& name)
    : event_out{f, std::move(ohepmc), std::move(ogemc), std::move(osimc),
                name} {
  create_branches();
}

void lA_out::push(const std::vector<lA_event>& events) {
  for (const auto& e : events) {
    push(e);
  }
}
void lA_out::push(const lA_event& e) {
  W_ = static_cast<float>(e.W());
  Q2_ = static_cast<float>(e.Q2());
  nu_ = static_cast<float>(e.nu());
  x_ = static_cast<float>(e.x());
  y_ = static_cast<float>(e.y());
  epsilon_ = static_cast<float>(e.epsilon());
  R_ = static_cast<float>(e.R());
  t_ = static_cast<float>(e.t());
  xv_ = static_cast<float>(e.xv());
  Q2plusM2_ = static_cast<float>(e.Q2plusM2());
  target_index_ = static_cast<int16_t>(e.target_index());
  photon_index_ = static_cast<int16_t>(e.photon_index());
  scat_index_ = static_cast<int16_t>(e.scat_index());
  leading_index_ = static_cast<int16_t>(e.leading_index());
  recoil_index_ = static_cast<int16_t>(e.recoil_index());

  rc_W_ = static_cast<float>(e.detected_W());
  rc_Q2_ = static_cast<float>(e.detected_Q2());
  rc_nu_ = static_cast<float>(e.detected_nu());
  rc_x_ = static_cast<float>(e.detected_x());
  rc_y_ = static_cast<float>(e.detected_y());
  rc_t_ = static_cast<float>(e.detected_t());
  rc_xv_ = static_cast<float>(e.detected_xv());
  rc_Q2plusM2_ = static_cast<float>(e.detected_Q2plusM2());
  rc_photon_index_ = static_cast<int16_t>(e.detected_photon_index());
  rc_scat_index_ = static_cast<int16_t>(e.detected_scat_index());
  rc_leading_index_ = static_cast<int16_t>(e.detected_leading_index());
  rc_recoil_index_ = static_cast<int16_t>(e.detected_recoil_index());

  event_out::push(e);
}

void lA_out::create_branches() {
  tree()->Branch("W", &W_);
  tree()->Branch("Q2", &Q2_);
  tree()->Branch("nu", &nu_);
  tree()->Branch("x", &x_);
  tree()->Branch("y", &y_);
  tree()->Branch("epsilon", &epsilon_);
  tree()->Branch("R", &R_);
  tree()->Branch("t", &t_);
  tree()->Branch("xv", &xv_);
  tree()->Branch("Q2plusM2", &Q2plusM2_);
  tree()->Branch("target_index", &target_index_);
  tree()->Branch("photon_index", &photon_index_);
  tree()->Branch("scat_index", &scat_index_);
  tree()->Branch("leading_index", &leading_index_);
  tree()->Branch("recoil_index", &recoil_index_);

  tree()->Branch("rc_W", &rc_W_);
  tree()->Branch("rc_Q2", &rc_Q2_);
  tree()->Branch("rc_nu", &rc_nu_);
  tree()->Branch("rc_x", &rc_x_);
  tree()->Branch("rc_y", &rc_y_);
  tree()->Branch("rc_t", &rc_t_);
  tree()->Branch("rc_xv", &rc_xv_);
  tree()->Branch("rc_Q2plusM2", &rc_Q2plusM2_);
  tree()->Branch("rc_photon_index", &rc_photon_index_);
  tree()->Branch("rc_scat_index", &rc_scat_index_);
  tree()->Branch("rc_leading_index", &rc_leading_index_);
  tree()->Branch("rc_recoil_index", &rc_recoil_index_);
}

} // namespace lager
