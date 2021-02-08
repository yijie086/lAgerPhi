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

#include "data.hh"
#include <cmath>
#include <lager/physics/photon.hh>

namespace lager {
namespace initial {

// =======================================================================================
// Create a new real target lager::particle
// =======================================================================================
target target::make_target(const lager::particle& ion, lager::particle nucl,
                           std::vector<lager::particle> remn, const double xs) {
  // get boost from ion frame to lab frame
  lager::particle::Boost boost_to_lab{-ion.p().BoostToCM()};
  // target_cm and remannt_cm are given in the ion helicity frame
  // ===> 1. rotate to the ion rest frame with the z-axis parallel to the lab
  // frame
  //      2. boost to lab frame
  nucl.rotate_uz(ion.p());
  nucl.boost(boost_to_lab);
  std::for_each(remn.begin(), remn.end(), [&](lager::particle& part) {
    part.rotate_uz(ion.p());
    part.boost(boost_to_lab);
    part.update_status(particle::status_code::SPECTATOR);
    part.vertex() = ion.vertex();
  });
  // now create and dress our target data
  target data{nucl.p(), nucl.type(), xs};
  data.particle().vertex() = ion.vertex();
  data.remnant_ = std::move(remn);

  return data;
}

// =======================================================================================
// Create a new real photon
// =======================================================================================
photon photon::make_real(const lager::particle& lepton,
                         const lager::particle& target, const double E,
                         const double xs = 1.) {
  lager::particle::XYZVector vec{lepton.p().Vect().Unit()};
  vec *= E;
  photon gamma{{vec.X(), vec.Y(), vec.Z(), E}, xs};
  gamma.scat_ = {lepton.type(), lepton.p() - gamma.particle().p(),
                 lager::particle::status_code::SCAT};
  gamma.W2_ = (gamma.particle().p() + target.p()).M2();
  gamma.nu_ = (gamma.particle().p()).Dot(target.p()) / target.mass();
  gamma.y_ = gamma.y_ * target.mass() / (lepton.p()).Dot(target.p());
  // this is not technically correct for a real photon, as the BS reaction might
  // have happened at any point earlier, but effectively this is correct near
  // the BS peak where the everything happens in the photon direction.
  gamma.particle().vertex() = lepton.vertex();
  gamma.scat_.vertex() = lepton.vertex();

  return gamma;
}

// =======================================================================================
// Create a new virtual photon
// =======================================================================================
photon photon::make_virtual(const lager::particle& lepton,
                            const lager::particle& target, const double Q2,
                            const double y, const double xs, const double phi) {

  // calculate scattered lepton
  // work in target rest frame
  lager::particle::Boost boost{target.p().BoostToCM()};
  lager::particle beam = lepton;
  beam.boost(boost);
  // now work with z-axis along the lepton beam direction
  const double E0 = beam.energy();
  const double E1 = E0 * (1. - y);
  const double theta1 = 2 * asin(sqrt(Q2 / (4. * E0 * E1)));
  const double phi1 = phi;
  lager::particle scat{
      lepton.type(),
      {cos(phi1) * sin(theta1), sin(phi1) * sin(theta1), cos(theta1)},
      E1,
      lager::particle::status_code::SCAT};
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
  vphoton.particle().vertex() = lepton.vertex();
  vphoton.scat_.vertex() = lepton.vertex();

  return vphoton;
}

} // namespace initial
} // namespace lager
