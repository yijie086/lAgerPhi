#include "photon.hh"
#include <cmath>
#include <liege/physics/photon.hh>

namespace liege {
namespace beam {

// =======================================================================================
// Create a new real photon
// =======================================================================================
photon photon::make_real(const particle& lepton, const particle& target,
                         const double E, const double xs = 1.) {
  particle::XYZVector vec{lepton.p().Vect().Unit()};
  vec *= E;
  photon gamma{{vec.X(), vec.Y(), vec.Z(), E}, xs};
  gamma.scat_ = {lepton.type(), lepton.p() - gamma.beam().p(),
                 particle::status_code::SCAT};
  gamma.W2_ = (gamma.beam().p() + target.p()).M2();
  gamma.nu_ = (gamma.beam().p()).Dot(target.p()) / target.mass();
  gamma.y_ = gamma.y_ * target.mass() / (lepton.p()).Dot(target.p());
  // this is not technically correct for a real photon, as the BS reaction might
  // have happened at any point earlier, but effectively this is correct near
  // the BS peak where the everything happens in the photon direction.
  gamma.beam().vertex() = lepton.vertex();
  gamma.scat_.vertex() = lepton.vertex();

  return gamma;
}

// =======================================================================================
// Create a new virtual photon
// =======================================================================================
photon photon::make_virtual(const particle& lepton, const particle& target,
                            const double Q2, const double y, const double xs,
                            const double phi) {

  // calculate scattered lepton
  // work in target rest frame
  particle::Boost boost{target.p().BoostToCM()};
  particle beam = lepton;
  beam.boost(boost);
  // now work with z-axis along the lepton beam direction
  const double E0 = beam.energy();
  const double E1 = E0 * (1. - y);
  const double theta1 = 2 * asin(sqrt(Q2 / (4. * E0 * E1)));
  const double phi1 = phi;
  particle scat{lepton.type(),
                {cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)},
                E1,
                particle::status_code::SCAT};
  // rotate to regular target rest frame, then boost to lab frame
  scat.rotate_uz(beam.p());
  scat.boost(boost.Inverse());

  // calculate the actual photon 4-momentum
  photon vphoton{lepton.p() - scat.p(), xs};
  vphoton.scat_ = scat;

  // epsilon
  vphoton.epsilon_ = physics::epsilon(Q2, y, lepton, target);

  // invariants
  vphoton.Q2_ = Q2;
  vphoton.y_ = y;
  vphoton.nu_ = y * (target.p()).Dot(lepton.p()) / target.mass();
  vphoton.W2_ = target.mass2() + 2 * target.mass() * vphoton.nu_ - Q2;
  vphoton.x_ = Q2 / (2 * target.mass() * vphoton.nu_);
  vphoton.beam().vertex() = lepton.vertex();
  vphoton.scat_.vertex() = lepton.vertex();

  return vphoton;
}

} // namespace beam
} // namespace liege
