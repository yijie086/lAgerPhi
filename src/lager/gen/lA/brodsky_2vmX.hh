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

#ifndef LAGER_GEN_LA_BRODSKY_2VMX_LOADED
#define LAGER_GEN_LA_BRODSKY_2VMX_LOADED

#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/lA/generator.hh>
#include <lager/gen/lA_event.hh>

namespace lager {
namespace lA {

// =============================================================================
// lA::brodsky_2vmX
//
// gamma + p -> VM + X process
//
// Uses the following expressions (cf. lager/physics/vm.hh)
//  * R (sigma_L/sigma_T):
//        R_vm_martynov(...)
//  * Dipole FF for sigma_gamma -> sigma_t:
//        dipole_ff_vm_hermes(...)
//  * t-channel cross section:
//        dsigma_dexp_bt_brodsky(...)
// =============================================================================
class brodsky_2vmX : public lA::generator {
public:
  using base_type = lA::generator;

  brodsky_2vmX(const configuration& cf, const string_path& path,
               std::shared_ptr<TRandom> r);
  virtual lA_event generate(const lA_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return max_exp_bt_range_.width(); }

private:
  double calc_max_xsec(const configuration& cf) const;
  interval<double> calc_max_t_range(const configuration& cf) const;

  // kinematics range
  interval<double> exp_bt_range(const double W2, const double Q2,
                                const double Mt) const;

  // cross section component evaluation
  double dsigma_dexp_bt(const double W2, const double Mt) const;
  double R(const double Q2) const;
  double dipole(const double Q2) const;

  // jacobian for d/dexp_bt -> d/dt
  double jacobian(const double t) const;

  // threshold squared for these particular particles (correctly handels the
  // case of particles with non-zero width)
  double threshold2(const particle& vm, const particle& recoil) const;

  // utility function
  lA_event make_event(const lA_data& initial, const double t, particle vm1,
                      particle X1, const double xs, const double R);

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;

  // cross section settings
  const double photo_b_;   // target FF constant
  const double photo_c2g_; // 2-gluon amplitude norm
  const double photo_c3g_; // 3-gluon amplitude norm
  const double R_vm_c_;    // c-parameter for R
  const double R_vm_n_;    // n-parameter for R
  const double dipole_n_;  // n-parameter for dipole factor

  // t-range and cross setion maxima
  const interval<double> max_t_range_;
  const interval<double> max_exp_bt_range_;
  const double max_;
};

} // namespace lA
} // namespace lager

#endif
