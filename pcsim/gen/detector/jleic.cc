#include "jleic.hh"

namespace pcsim {
namespace detect {
FACTORY_REGISTER(detector, jleic, "jleic");

bool jleic::accepted(const particle& part) {

  const double kin[4] = {part.momentum(), part.p().theta() * TMath::RadToDeg(),
                         part.p().phi() * TMath::RadToDeg()};
  int region[1];
  double acc[1];
  acceptance_.get_acceptance(part.type<int>(), kin, region, acc, false);
  return (rng()->Uniform(0, 1) < acc[0]);
}
detected_particle jleic::generate(const particle& part) {
  const double kin[4] = {part.momentum(), part.p().theta() * TMath::RadToDeg(),
                         part.p().phi() * TMath::RadToDeg()};
  double res[4];
  double smeared[4];
  int region[1];
  double acc[1];
  acceptance_.get_acceptance(part.type<int>(), kin, region, acc, false);
  resolution_.get_resolution(part.type<int>(), kin, region[0], res, smeared,
                             false);
  particle::Polar3DVector p3s{smeared[0], smeared[1] * TMath::DegToRad(),
                              smeared[2] * TMath::DegToRad()};
  return {part, {p3s.X(), p3s.Y(), p3s.Z(), smeared[3]}, region[0]};
}

} // namespace detect
} // namespace pcsim
