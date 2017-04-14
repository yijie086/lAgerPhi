#include "decay.hh"
#include <cmath>
#include <pcsim/core/pdg.hh>

namespace pcsim {
namespace decay {

void decay_handler::process(event& e, const double epsilon_R) {
  for (int i = 0; i < e.size(); ++i) {
    // we won't decay stable particles
    if (e[i].stable()) {
      continue;
    }
    // SCHC leptonic e+e- decay of vms
    if (e[i].type() == pdg_id::J_psi || e[i].type() == pdg_id::upsilon ||
        e[i].type() == pdg_id::phi) {
      e.update_weight(e[i].pdg()->DecayChannel(0)->BranchingRatio());
      std::pair<particle, particle> decay_products{{pdg_id::e_plus},
                                                   {pdg_id::e_minus}};
      const double r04 = epsilon_R / (1 + epsilon_R);
      const double phi = rng()->Uniform(0., TMath::TwoPi());
      const double ctheta =
          rand_f({-1, 1},
                 [=](const double ctheta) {
                   return ((1. + r04) + (1. - 3. * r04) * ctheta * ctheta);
                 },
                 2.001);
      const double theta = acos(ctheta);
      decay_2body(e[i], theta, phi, decay_products);
      e.add_daughter(decay_products.first, i);
      e.add_daughter(decay_products.second, i);
      e[i].update_status(particle::status_code::DECAYED_SCHC);
    }
    // Pc decay according to wang et al.
    else if (e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 5, 1, 1) ||
             e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 5, 1, 4) ||
             e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 5, -1, 2) ||
             e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 5, -1, 3) ||
             e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 3, 1, 3) ||
             e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 3, 1, 2) ||
             e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 3, -1, 4) ||
             e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 3, -1, 1)) {
      std::pair<particle, particle> decay_products{
          {pdg_id::J_psi, particle::status_code::UNSTABLE_SCHC}, {pdg_id::p}};
      const double phi = rng()->Uniform(0., TMath::TwoPi());
      double ctheta = -1;
      if (e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 5, 1, 1) ||
          e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 5, 1, 4)) {
        // result from a pol6 fit to a digitized version of figure 6c from
        // PRD92-034022(2015)
        ctheta = rand_f({-1, 1},
                        [](const double x) {
                          const double x2 = x * x;
                          const double x3 = x2 * x;
                          const double x4 = x3 * x;
                          const double x5 = x4 * x;
                          const double x6 = x5 * x;
                          return .149211 - 0.194418 * x - 0.563191 * x2 +
                                 0.374024 * x3 + 0.658942 * x4 + 0.110057 * x5 +
                                 0.0931712 * x6;
                        },
                        0.63);
      } else if (e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 5, -1, 2) ||
                 e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 5, -1, 3)) {
        // result from a pol7 fit to a digitized version of figure 5c from
        // PRD92-034022(2015)
        ctheta = rand_f({-1, 1},
                        [](const double x) {
                          const double x2 = x * x;
                          const double x3 = x2 * x;
                          const double x4 = x3 * x;
                          const double x5 = x4 * x;
                          const double x6 = x5 * x;
                          const double x7 = x6 * x;
                          return 1.31241 - 1.19802 * x + 1.58351 * x2 +
                                 17.1514 * x3 + 20.8306 * x4 - 4.43848 * x5 +
                                 2.67151 * x6 + 6.06378 * x7;
                        },
                        44.06);
      } else if (e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 3, 1, 3) ||
                 e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 3, 1, 2)) {
        // result from a expo fit to a digitized version of figure 5b from
        // PRD92-034022(2015)
        ctheta = rand_f({-1, 1}, [](const double x) { return exp(-5.944 - x); },
                        0.00713);
      } else if (e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 3, -1, 4) ||
                 e[i].type() == pdg_pentaquark_id(4, 2, 2, 1, 4, 3, -1, 1)) {
        // result from a pol2 fit to a digitized version of figure 6b from
        // PRD92-034022(2015)
        ctheta = rand_f({-1, 1},
                        [](const double x) {
                          const double x2 = x * x;
                          return 0.00845846 - 0.0128146 * x + 0.00526053 * x2;
                        },
                        0.0266);
      }
      const double theta = acos(ctheta);
      decay_2body(e[i], theta, phi, decay_products);
      e.add_daughter(decay_products.first, i);
      e.add_daughter(decay_products.second, i);
      e[i].update_status(particle::status_code::DECAYED);
    }
  }
  // that's all
}

// two body decay of a particle 'part' into two particles xx (xx.first,
// xx.second), with angles of the first decay particle ('theta1', 'phi1')
//
// Note: the angles are assumed to be in the helicity frame of 'part'
void decay_2body(const particle& part, const double theta_1, const double phi_1,
                 std::pair<particle, particle>& xx) {

  // calculate the CM kinematics
  const double E = part.mass();
  const double M2_1 = xx.first.mass() * xx.first.mass();
  const double M2_2 = xx.second.mass() * xx.second.mass();
  const double E_1 = (E * E + M2_1 - M2_2) / (2 * E);
  const double E_2 = (E * E - M2_1 + M2_2) / (2 * E);
  const double P_1 = sqrt(E_1 * E_1 - M2_1);
  const double P_2 = sqrt(E_2 * E_2 - M2_2);

  // create the 4-vectors in the part helicity frame
  TVector3 mom;
  mom.SetMagThetaPhi(P_1, theta_1, phi_1);
  xx.first.p() = {mom, E_1};
  // decay particle 1 and two are back-to-back in the CM frame
  mom.SetMagThetaPhi(P_2, theta_1, phi_1);
  xx.second.p() = {-mom, E_2};

  // boost back to the rotated version of the original frame pointing in the
  // direction of part
  const particle::Boost b{
      -(particle::XYZTVector{0, 0, part.p().Vect().Mag(), part.p().E()}
            .BoostToCM())};
  xx.first.boost(b);
  xx.first.boost(b);

  // rotate to the original frame
  xx.first.rotate_uz(part.p().Vect());
  xx.second.rotate_uz(part.p().Vect());

  // that's all
}

} // decay
} // pcsim
