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

#ifndef LAGER_GEN_INITIAL_DATA_LOADED
#define LAGER_GEN_INITIAL_DATA_LOADED

#include <lager/core/generator.hh>
#include <lager/core/particle.hh>

namespace lager {
namespace initial {

// VERTEX DATA
using vertex = particle::XYZTVector;

// PRIMARY BEAM DATA
// =============================================================================
// initial::beam
//
// primary beam data
// =============================================================================
class beam : public generator_data {
public:
  beam() = default;
  beam(const particle& part) : part_{part} {}
  beam(const particle& part, const double xs)
      : generator_data{xs}, part_{part} {}
  beam(const double xs) : generator_data{xs} {}

  lager::particle& particle() { return part_; }
  const lager::particle& particle() const { return part_; }

private:
  lager::particle part_;
};

// SECONDARY BEAM DATA
// =============================================================================
// initial::target
//
// secondary target beam in a primary ion beam (e.g. beam proton,
// or proton in nucleon)
// =============================================================================
class target : public beam {
public:
  // generalized constructor:
  // split a nucleon beam in a nucleon jet and a remnant
  //
  // nucl and remn are given in the ion helicity frame
  static target make_target(const lager::particle& ion, lager::particle nucl,
                            std::vector<lager::particle> remn = {},
                            const double xs = 1);

  target() = default;
  target(const target&) = default;
  target& operator=(const target&) = default;

  target(const double xs) : beam{xs} {}

  target(const particle::XYZTVector& p, const pdg_id& id)
      : beam{{id, p, lager::particle::status_code::SECONDARY_BEAM}} {}
  target(const particle::XYZTVector& p, const pdg_id& id, const double xs)
      : beam{{id, p, lager::particle::status_code::SECONDARY_BEAM}, xs} {}

  const std::vector<lager::particle>& remnant() const { return remnant_; }

private:
  std::vector<lager::particle> remnant_; // nucleon remnant (if any)
};

// =============================================================================
// initial::photon
//
// secondary photon beam (e.g. virtual photon, or bremsstrahlung photon)
// =============================================================================
class photon : public beam {
public:
  // generalized constructors:
  // make collinear real photon event with energy E
  static photon make_real(const lager::particle& lepton,
                          const lager::particle& target, const double E,
                          const double xs);

  // generate virtual photon event with Q2 and y
  // also needs an azimuthal angle in the CM frame
  static photon make_virtual(const lager::particle& lepton,
                             const lager::particle& target, const double Q2,
                             const double y, const double xs, const double phi);

  photon() = default;
  photon(const photon&) = default;
  photon& operator=(const photon&) = default;

  photon(const double xs) : beam{xs} {}

  photon(const lager::particle::XYZTVector& p)
      : beam{{pdg_id::gamma, p, lager::particle::status_code::SECONDARY_BEAM}} {
  }
  photon(const lager::particle::XYZTVector& p, const double xs)
      : beam{{pdg_id::gamma, p, lager::particle::status_code::SECONDARY_BEAM},
             xs} {}

  double epsilon() const { return epsilon_; }
  double W() const { return sqrt(W2_); }
  double W2() const { return W2_; }
  double Q2() const { return Q2_; }
  double nu() const { return nu_; }
  double x() const { return x_; }
  double y() const { return y_; }

  lager::particle& scat() { return scat_; }
  const lager::particle& scat() const { return scat_; }

private:
  double epsilon_{0.};   // gamma_L/gamma_T (xs stores just gamma_T)
  double W2_{0.};        // invariant mass of photon-target system
  double Q2_{0.};        // photon virtuality
  double nu_{0.};        // photon enery in target rest frame
  double x_{0.};         // Bjorken x
  double y_{0.};         // energy fraction of photon in target rest frame
  lager::particle scat_; // scattered lepton
};
} // namespace initial
} // namespace lager

#endif
