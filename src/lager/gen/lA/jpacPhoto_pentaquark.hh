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

#ifndef LAGER_GEN_LA_JPACPHOTO_PENTAQUARK_LOADED
#define LAGER_GEN_LA_JPACPHOTO_PENTAQUARK_LOADED

#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/lA/generator.hh>
#include <lager/gen/lA_event.hh>

#include <jpacPhoto/amplitudes/amplitude_sum.hpp>
#include <jpacPhoto/amplitudes/baryon_resonance.hpp>
#include <jpacPhoto/reaction_kinematics.hpp>

namespace lager {
namespace lA {

// =============================================================================
// lA::jpacPhoto_pentaquark
//
// gamma + p -> VM + p process
// (using pentaquark exchange)
//
// (uses additional dipole formfactor to model Q2 dependence)
//
// =============================================================================
class jpacPhoto_pentaquark : public lA::generator {
public:
  using base_type = lA::generator;

  jpacPhoto_pentaquark(const configuration& cf, const string_path& path,
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
  std::vector<std::unique_ptr<jpacPhoto::baryon_resonance>>
  init_components() const;
  jpacPhoto::amplitude_sum init_ampl();

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;
  const particle decay_lplus_;
  const particle decay_lminus_;

  // cross section settings
  const std::vector<int> spin_parity_;  // resonance (spin x 2) times parity
  const std::vector<double> mass_;      // mass of resonance
  const std::vector<double> width_;     // width of resonance
  const std::vector<std::string> name_; // resonance names
  const double photo_balance_;          // balance between photocouplings
  const double branching_;              // branching fraction to J/psi-p
  const double
      dipole_n_; // n-parameter for dipole factor for slightly improved EPA

  // cross section utility variables
  std::unique_ptr<jpacPhoto::reaction_kinematics>
      reaction_; // Reaction kinematics
  std::vector<std::unique_ptr<jpacPhoto::baryon_resonance>> components_;
  jpacPhoto::amplitude_sum ampl_;

  // t-range and cross section maxima
  const interval<double> max_t_range_;
  const double max_;
};

} // namespace lA
} // namespace lager

#endif
