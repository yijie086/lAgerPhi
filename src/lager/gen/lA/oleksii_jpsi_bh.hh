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

#ifndef LAGER_GEN_LA_OLEKSII_JPSI_BH_LOADED
#define LAGER_GEN_LA_OLEKSII_JPSI_BH_LOADED

#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/lA/generator.hh>
#include <lager/gen/lA_event.hh>

namespace lager {
namespace lA {

// =============================================================================
// lA::oleksii_jpsi_bh
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
class oleksii_jpsi_bh : public lA::generator {
public:
  using base_type = lA::generator;

  oleksii_jpsi_bh(const configuration& cf, const string_path& path,
                  std::shared_ptr<TRandom> r);
  virtual lA_event generate(const lA_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const {
    // last factor is the cos(theta) range
    // return max_t_range_.width() *
    return (exp(1.13 * max_t_range_.max) - exp(1.13 * max_t_range_.min)) *
           (Mll_range_.max * Mll_range_.max - Mll_range_.min * Mll_range_.min) *
           TMath::TwoPi() * 2.;
  }

private:
  interval<double> calc_max_t_range(const configuration& cf) const;

  // threshold squared for these particular particles (correctly handels the
  // case of particles with non-zero width)
  double threshold2(const particle& vm, const particle& recoil) const;

  // utility function
  lA_event make_event(const lA_data& initial, const double t,
                            particle vm1, particle X1, const double xs,
                            const double thetaCM_el,
                            const double phiCM_el) const;

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;

  // t-range and cross setion maxima
  const interval<double> max_t_range_;
  const interval<double> Mll_range_;
  const double max_;

  // subtraction constant
  const double T_0_;

  // theta acceptance to cut out colinear enhancements
  const interval<double> p_range_;
  const interval<double> theta_range_;
};

} // namespace lA
} // namespace lager

#endif
