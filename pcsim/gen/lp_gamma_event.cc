#include "lp_gamma_event.hh"

namespace pcsim {

// =============================================================================
// IMPLEMENTATION: event_out
// =============================================================================
lp_gamma_out::lp_gamma_out(std::shared_ptr<TFile> f, const std::string& name)
    : event_out{f, name} {
  create_branches();
}

void lp_gamma_out::push(const std::vector<lp_gamma_event>& events) {
  for (const auto& e : events) {
    push(e);
  }
}
void lp_gamma_out::push(const lp_gamma_event& e) {
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

void lp_gamma_out::create_branches() {
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

} // namespace pcsim
