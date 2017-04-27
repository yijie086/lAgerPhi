#define JLEIC_IMPL
#include "jleic.hh"

namespace pcsim {
namespace detector {

void jleic::process(event& e) const {
  for (auto& part : e) {
    if (part.final_state()) {
      const double kin[4] = {part.momentum(),
                             part.p().theta() * TMath::RadToDeg(),
                             part.p().phi() * TMath::RadToDeg()};
      int region[1];
      double acc[1];
      acceptance_.get_acceptance(part.type<int>(), kin, region, acc, false);
      if (rng()->Uniform(0, 1) > acc[0]) {
        // not accepted, continue
        continue;
      }
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
