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

#ifndef LAGER_GEN_LA_HOLOGRAPHIC_VM_LOADED
#define LAGER_GEN_LA_HOLOGRAPHIC_VM_LOADED

#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/lA/generator.hh>
#include <lager/gen/lA_event.hh>

namespace lager {
namespace lA {

// =============================================================================
// lA::holographic_vm
//
// gamma/gamma* + p -> VM + X process
//
// Uses the following expressions (cf. lager/physics/vm.hh)
//  * dsigma_dt_holographic --> photoproduction cross section
//  * Dipole Q2 dependence for sigma_gamma -> sigma_t:
//        dipole_ff_vm_hermes(...)
//  * R (sigma_L/sigma_T):
//        R_vm_martynov(...)
//
// =============================================================================
class holographic_vm : public lA::generator {
public:
  using base_type = lA::generator;

  holographic_vm(const configuration& cf, const string_path& path,
                 std::shared_ptr<TRandom> r);
  virtual lA_event generate(const lA_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return max_t_range_.width(); }

private:
  double calc_max_xsec(const configuration& cf) const;
  interval<double> calc_max_t_range(const configuration& cf,
                                    const string_path& path) const;

  // FIXME this should be a general utility function
  // threshold squared for these particular particles (correctly handels the
  // case of particles with non-zero width)
  double threshold2(const particle& vm, const particle& recoil) const;

  // utility function
  lA_event make_event(const lA_data& initial, const double t, particle vm1,
                      particle X1, const double xs, const double R);

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;

  // total cross section settings
  const double A0_;
  const double m_A_;
  const double C0_;
  const double m_C_;
  const double N_;
  // R settings
  const double R_vm_c_; // c-parameter for R
  const double R_vm_n_; // n-parameter for R
  // Q2 dependence settings
  const double dipole_n_; // n-parameter for dipole factor

  // t-range and cross setion maxima
  const interval<double> max_t_range_;
  const double max_;
};

} // namespace lA
} // namespace lager

#endif
