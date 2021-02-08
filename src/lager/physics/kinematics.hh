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

#ifndef LAGER_PHYSICS_KINEMATICS_LOADED
#define LAGER_PHYSICS_KINEMATICS_LOADED

#include <cmath>
#include <lager/core/interval.hh>
#include <lager/core/particle.hh>

namespace lager {
namespace physics {

// =======================================================================================
// calculate the allowed Q2 range for a given lepton beam, target and y
//
// The lower bound is given by the kinematics of t-channel scattering,
// while the upper bound is given by either t-channel kinematics, or the
// requirement that the final W be larger than the target mass
// =======================================================================================
inline interval<double> Q2_range(const particle& beam, const particle& target,
                                 const double y) {
  const double E = (beam.p()).Dot(target.p()) / target.mass();
  const double E2 = E * E;
  const double m2 = beam.mass2();
  // lower and upper bound from t-channel process on electron leg
  const double comp1 = -2. * m2 + 2. * E2 * (1. - y);
  const double comp2 = 2. * sqrt((E2 - m2) * (E2 * (1. - y) * (1. - y) - m2));
  const double Q2_low = comp1 - comp2;
  const double Q2_high1 = comp1 + comp2;
  // alternative upper bound from requirement that final state has at least the
  // invariant mass of the target mass (W2min = target.mass^2, meaning that 2 M
  // nu = Q2)
  const double Q2_high2 = 2 * target.mass() * E * y;
  return {Q2_low, fmin(Q2_high1, Q2_high2)};
}

// =======================================================================================
// calculate the allowed t range
//
// input:
//   * W2 (=s)
//   * Q2 (-photon mass squared)
//   * Mt (target mass)
//   * Mv (VM mass)
//   * Mr (recoil mass)
// =======================================================================================
inline interval<double> t_range(const double W2, const double Q2,
                                const double Mt, const double Mv,
                                const double Mr) {
  const double Mt2 = Mt * Mt;
  const double Mr2 = Mr * Mr;
  const double Mv2 = Mv * Mv;
  const double W = sqrt(W2);
  // CM relations for energy and momenta
  const double Et_cm = (W2 + Q2 + Mt2) / (2. * W);
  const double Pt_cm = sqrt(Et_cm * Et_cm - Mt2);
  const double Er_cm = (W2 - Mv2 + Mr2) / (2. * W);
  const double Pr_cm = sqrt(Er_cm * Er_cm - Mr2);
  // t_low when recoil direction changes by 90 degrees compared to target,
  // t_high when the recoil in the target direction
  // note: what we call t_high is also refered to as t_min
  const double t_low = Mt2 + Mr2 - 2. * Et_cm * Er_cm - 2. * Pt_cm * Pr_cm;
  const double t_high = Mt2 + Mr2 - 2. * Et_cm * Er_cm + 2. * Pt_cm * Pr_cm;
  // that's all!
  return {t_low, t_high};
}

} // namespace physics
} // namespace lager

#endif
