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

#ifndef LAGER_GEN_INITIAL_TARGET_GEN_LOADED
#define LAGER_GEN_INITIAL_TARGET_GEN_LOADED

#include <TRandom.h>
#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/initial/data.hh>
#include <lager/gen/initial/generator.hh>
#include <memory>

namespace lager {
namespace initial {

// estimated target info from config file to be used by lA generators to
// estimate phase-space and cross section ranges
particle estimated_target(const configuration& cf);

// simple target identical to the ion beam
class primary_target : public target_generator {
public:
  primary_target(const configuration& cf, const string_path& path,
                 std::shared_ptr<TRandom> r);

  virtual target generate(const beam& ion);
  // no contribution to the cross section and phase space
  virtual double max_cross_section() const { return 1; }
  virtual double phase_space() const { return 1; }

private:
  const particle target_;
};

// Fermi distribution using old fortran routines from J.S O'Connell and J.W.
// Lightbody Jr. See physics/fermi.hh for more info
// Note that this routine is very generic, and we do not treat remnant particles
// here
class fermi87 : public target_generator {
public:
  fermi87(const configuration& cf, const string_path& path,
          std::shared_ptr<TRandom> r);

  virtual target generate(const beam& ion);
  // We treat the nucleon selection as orthogonal to the main generation
  // process, so we first select a nucleon according to the fermi-distribution.
  // Because of this we do not have to cary the cross section and phase-space
  // information into the main accept-reject procedure.
  // I expect this to be more computationally efficient in most cases (although
  // I did not
  // test the difference).
  virtual double max_cross_section() const { return 1; }
  virtual double phase_space() const { return 1; }

private:
  // calculate the normalization to the fermi-distribution function
  double calc_norm() const;
  // calculate the maximum fermi-distribution PDF
  double calc_max() const;
  // Actual p.d.f implementation
  double pdf(const double P) const;

  const int A_;            // nucleus A
  const particle nucleon_; // nucleon-of-interest
  const double k_max_;     // max fermi momentum
  const double norm_;      // internal normalization factor
  const double xs_max_;    // maximum value of the PDF
};

} // namespace initial
} // namespace lager

#endif
