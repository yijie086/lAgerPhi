#define JLEIC_IMPL
#include "jleic.hh"

namespace pcsim {
namespace detector {

void jleic::process(event& e) const {
  for (auto& part : e) {
    if (part.final_state()) {
      LOG_JUNK2("jleic", "Checking acceptance for final state particle " +
                             part.name() + " (status: " +
                             std::to_string(part.status<int>()) + ")");
      const double kin[4] = {part.momentum(),
                             part.p().theta() * TMath::RadToDeg(),
                             part.p().phi() * TMath::RadToDeg(), part.energy()};
      int region[1];
      double acc[1];
      acceptance_.get_acceptance(part.type<int>(), kin, region, acc, false);
      LOG_JUNK2("jleic",
                "Calculated acceptance: " + std::to_string(acc[0] * 100) +
                    "\% (region: " + std::to_string(region[0]) + ")");
      if (rng()->Uniform(0, 1) > acc[0]) {
        LOG_JUNK2("jleic", "Not accepted");
        // not accepted, continue
        continue;
      }
      LOG_JUNK2("jleic", "Particle accepted!");
      double res[4];
      double smeared[4];
      resolution_.get_resolution(part.type<int>(), kin, region[0], res, smeared,
                                 false);
      particle::Polar3DVector p3s{smeared[0], smeared[1] * TMath::DegToRad(),
                                  smeared[2] * TMath::DegToRad()};
      e.add_detected(
          {part, {p3s.X(), p3s.Y(), p3s.Z(), smeared[3]}, region[0]});
    }
  }
}

} // namespace detector
} // namespace pcsim
