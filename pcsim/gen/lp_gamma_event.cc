#include "lp_gamma_event.hh"

namespace pcsim {

// =============================================================================
// IMPLEMENTATION: event_out
// =============================================================================
lp_gamma_out::lp_gamma_out(std::shared_ptr<TFile> f, const std::string& name)
    : event_out{f, name} {
  create_branches();
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
  tree()->Branch("photon", &photon_index_);
  tree()->Branch("scat", &scat_index_);
  tree()->Branch("leading", &leading_index_);
  tree()->Branch("recoil", &recoil_index_);
}

} // namespace pcsim
