#include "lp_gamma.hh"
#include <cmath>
#include <pcsim/core/particle.hh>
#include <pcsim/core/pdg.hh>

namespace pcsim {
namespace reconstruction {

void lp_gamma::process(lp_gamma_event& e) const {
  // do additional event reconstruction (reconstruction of scattered
  // lepton/photon/recoil kinematics are done automatically be
  // lp_gamma_event::add_detected())
  for (const auto& i = 0; i < e.detected().size(); ++i) {
    // leading particle from decay products
    if (e.detected_leading_index() < 0 && e.leading_index() >= 0) {
      if (e.detected(i).parent_first() == e.leading_index()) {
        // locate the other decay product
        for (int j = i + 1; j < e.detected().size(); ++j) {
          if (e.detected(j).parent_first() == e.leading_index()) {
            // we found both!
            e.add_detected(e.leading(), e.detected(i).p() + e.detected(j).p(),
                           -1);
          }
        }
      }
    }
  }
  // that's all
}

} // namespace reconstruction
} // namespace pcsim
