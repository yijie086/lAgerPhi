// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2020 Sylvester Joosten <sjoosten@anl.gov>
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

#ifndef LAGER_GEN_LA_OLEKSII_2VMp_LOADED
#define LAGER_GEN_LA_OLEKSII_2VMp_LOADED

#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/lA/generator.hh>
#include <lager/gen/lA_event.hh>

namespace lager {
namespace lA {

// =============================================================================
// lA::oleksii_2vmp
//
// gamma + p -> VM + p process
//
// Uses the following expressions (cf. lager/physics/vm.hh)
//  * R (sigma_L/sigma_T):
//        R_vm_martynov(...)
//  * Dipole FF for sigma_gamma -> sigma_t:
//        dipole_ff_vm_hermes(...)
//  * t-channel cross section:
//        Oleksii's implementation for J/psi and Upsilon
// =============================================================================

class oleksii_2vmp_amplitude;
class oleksii_2vmp_slope;

class oleksii_2vmp : public lA::generator {
public:
  using base_type = lA::generator;

  oleksii_2vmp(const configuration& cf, const string_path& path,
               std::shared_ptr<TRandom> r);
  virtual lA_event generate(const lA_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return max_exp_b0t_range_.width(); }

private:
  double calc_min_b() const;
  double calc_max_b(const configuration& cf) const;
  interval<double> calc_max_b_range(const configuration& cf) const;
  double calc_max_xsec(const configuration& cf) const;
  interval<double> calc_max_t_range(const configuration& cf) const;

  // kinematics range
  interval<double> t_range(const double W2, const double Q2,
                           const double Mt) const;
  // kinematics range
  interval<double> exp_bt_range(const double W2, const double Q2,
                                const double Mt) const;

  // cross section component evaluation
  double dsigma_dexp_b0t(const double W2, const double t, const double b) const;
  double R(const double Q2) const;
  double dipole(const double Q2) const;

  // jacobian (equal to unity)
  double jacobian(const double t) const;

  // threshold squared for these particular particles (correctly handels the
  // case of particles with non-zero width)
  double threshold2(const particle& vm, const particle& recoil) const;

  // utility function
  lA_event make_event(const lA_data& initial, const double t,
                            const double b, particle vm1, particle X1,
                            const double xs, const double R);

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;

  // scattering amplitude
  class oleksii_2vmp_amplitude* ampl_;
  class oleksii_2vmp_slope* slope_;

  // cross section settings
  const double T0_;       // subtraction constant for dispersion relation
  const double R_vm_c_;   // c-parameter for R
  const double R_vm_n_;   // n-parameter for R
  const double dipole_n_; // n-parameter for dipole factor

  // t-range and cross setion maxima
  const interval<double> max_b_range_; // upper limit to b parameter
  const interval<double> max_t_range_;
  const interval<double> max_exp_b0t_range_;
  const double max_;
};

} // namespace lA
} // namespace lager

#endif
