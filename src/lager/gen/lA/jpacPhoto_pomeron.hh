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

#ifndef LAGER_GEN_LA_JPACPHOTO_POMERON_LOADED
#define LAGER_GEN_LA_JPACPHOTO_POMERON_LOADED

#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/lA/generator.hh>
#include <lager/gen/lA_event.hh>

#include <jpacPhoto/amplitudes/pomeron_exchange.hpp>
#include <jpacPhoto/reaction_kinematics.hpp>

namespace lager {
namespace lA {

// =============================================================================
// lA::jpacPhoto_pomeron
//
// gamma + p -> VM + p process
// (using pomeron exchange)
//
// (uses additional dipole formfactor to model Q2 dependence)
//
// =============================================================================
class jpacPhoto_pomeron : public lA::generator {
public:
  using base_type = lA::generator;

  jpacPhoto_pomeron(const configuration& cf, const string_path& path,
                    std::shared_ptr<TRandom> r);
  virtual lA_event generate(const lA_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return max_t_range_.width(); }

private:
  double calc_max_xsec(const configuration& cf) /*const*/;
  interval<double> calc_max_t_range(const configuration& cf) const;

  // kinematics range
  interval<double> t_range(const double W2, const double Q2,
                           const double Mt) const;

  // cross section component evaluation
  // note:: non-const as we call a jpacPhoto function that's non-const
  double dsigma_dt(const double s, const double t) /*const*/
      ;
  double dipole(const double Q2) const;

  // threshold squared for these particular particles (correctly handels the
  // case of particles with non-zero width)
  double threshold2(const particle& vm, const particle& recoil) const;

  // utility function
  lA_event make_event(const lA_data& initial, const double t, particle vm1,
                      particle X1, const double xs);

  // initialize the reaction kinematics
  std::unique_ptr<jpacPhoto::reaction_kinematics> init_reaction() const;
  jpacPhoto::pomeron_exchange init_ampl();

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;
  const particle decay_lplus_;
  const particle decay_lminus_;

  // cross section settings
  const double regge_inter_; // Simple Regge trajectory intercept
  const double regge_slope_; // Simple Regge trajectory slope
  const double photo_norm_;  // t-channel normalization
  const double photo_slope_; // t-channel slope
  const double
      dipole_n_; // n-parameter for dipole factor for slightly improved EPA

  // cross section utility variables
  std::unique_ptr<jpacPhoto::reaction_kinematics>
      reaction_;                     // Reaction kinematics
  linear_trajectory regge_;          // Regge trajectory
  jpacPhoto::pomeron_exchange ampl_; // Amplitude

  // t-range and cross section maxima
  const interval<double> max_t_range_;
  const double max_;
};

} // namespace lA
} // namespace lager

#endif
