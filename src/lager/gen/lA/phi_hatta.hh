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

#ifndef LAGER_GEN_LA_PHI_HATTA_LOADED
#define LAGER_GEN_LA_PHI_HATTA_LOADED

#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/lA/generator.hh>
#include <lager/gen/lA_event.hh>

namespace lager {
namespace lA {

// =============================================================================
// lA::phi_hatta
//
// gamma/gamma* + p -> VM + X process
//
// Uses the following expressions (cf. lager/physics/vm.hh)
//  * sigma_T, including dipole-like Q2 dependence:
//        sigmaT_phi_clas
//  * (sigma_L/sigma_T):
//        R_phi_clas
//  * t-dependent dependence (normalized):
//        exp_ff_normalized or dipole_ff_normalized (configurable)
// =============================================================================
class phi_hatta : public lA::generator {
public:
  using base_type = lA::generator;

  phi_hatta(const configuration& cf, const string_path& path,
             std::shared_ptr<TRandom> r);
  virtual lA_event generate(const lA_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return max_t_range_.width(); }

private:
  double calc_max_xsec(const configuration& cf) const;
  interval<double> calc_max_t_range(const configuration& cf, const string_path& path) const;

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
  const double alpha_1_;
  const double alpha_2_;
  const double alpha_3_;
  const double nu_T_;
  // R settings
  const double c_R_;
  // FF function
  std::function<double(double /*Q2*/, double /*W*/, double /*t*/,
                       double /*Mt*/)>
      ff_func_;

  // t-range and cross setion maxima
  const interval<double> max_t_range_;
  const double max_;
};

} // namespace lA
} // namespace lager

#endif
