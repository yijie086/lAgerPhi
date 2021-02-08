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

#ifndef LAGER_GEN_LA_RESONANCE_QPQ_LOADED
#define LAGER_GEN_LA_RESONANCE_QPQ_LOADED

#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/lA/generator.hh>
#include <lager/gen/lA_event.hh>

namespace lager {
namespace lA {

// =============================================================================
// lA::resonance_qpq
//
// gamma + p -> quarkonium pentaquark resonance process
//
// Uses the following expressions (cf. lager/physics/vm.hh)
// (assuming the photo-production happens through a given quarkonium pole)
//  * R (sigma_L/sigma_T):
//        R_martynov(...)
//  * Dipole FF for sigma_gamma -> sigma_t:
//        dipole_ff__hermes(...)
//  * photo-production cross section:
//        simple resonance (BW)
// =============================================================================
class resonance_qpq : public lA::generator {
public:
  using base_type = lA::generator;
  resonance_qpq(const configuration& cf, const string_path& path,
                std::shared_ptr<TRandom> r);
  virtual lA_event generate(const lA_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return 1.; }

private:
  double calc_max_xsec(const configuration& cf) const;
  interval<double> calc_W2_range(const double n_sigma) const;
  double calc_normalization() const;

  // cross section component evaluation
  double sigma(const double W2) const;
  double R(const double Q2) const;
  double dipole(const double Q2) const;

  // recoil and  particle info
  const particle vm_pole_;          // relevant VM pole
  const particle qpq_;              // quarkonium pentaquark assumption

  // cross section settings
  const double mass_;      // resonance mass in GeV
  const double width_;     // Width of the resonance
  const double amplitude_; // peak cross section amplitude
  const double norm_;      // normalization for the BW (to have peak amplitude_)
  const double coupling_;  // coupling through the quarkonium pole
  const double R_vm_c_;    // c-parameter for R
  const double R_vm_n_;    // n-parameter for R
  const double dipole_n_;  // n-parameter for dipole factor

  const interval<double> W2_range_; // mass squared range around the pole mass


  const double max_;
};

} // namespace lA
} // namespace lager

#endif
