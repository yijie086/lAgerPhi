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

#ifndef LAGER_GEN_INITIAL_PHOTON_GEN_LOADED
#define LAGER_GEN_INITIAL_PHOTON_GEN_LOADED

#include <TRandom.h>
#include <lager/core/generator.hh>
#include <lager/core/particle.hh>
#include <lager/gen/initial/data.hh>
#include <lager/gen/initial/fixed_target.hh>
#include <lager/gen/initial/generator.hh>
#include <lager/physics/photon.hh>
#include <memory>

namespace lager {
namespace initial {

// no_photon generator, useful for e.g. DVCS + BH
class no_photon : public photon_generator {
public:
  no_photon(const configuration&, const string_path&, std::shared_ptr<TRandom>);

  virtual photon generate(const beam&, const target&);
  virtual double max_cross_section() const { return 1.; };
  virtual double phase_space() const { return 1.; }
};

// Bremsstrahlung photons
class bremsstrahlung : public photon_generator {
public:
  enum class model { FLAT, PARAM, APPROX, EXACT };

  bremsstrahlung(const configuration& cf, const string_path& path,
                 std::shared_ptr<TRandom> r);

  virtual photon generate(const beam&, const target&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return E_range_.width(); }

protected:
  double intensity(const double E, const double beam) const;

private:
  const model model_;   // bremsstrahlung model
  const double rl_;     // number of radiation lenghts (when using approx model,
                        // set to zero otherwise)
  const double E_beam_; // (maximum) electron beam energy
  const interval<double> E_range_; // photon energy range
  const double max_;               // the maximum intensity.
};

// Bremsstrahlung photons for a realistic (extended) target
class bremsstrahlung_realistic_target : public photon_generator {
public:
  bremsstrahlung_realistic_target(const configuration& cf,
                                  const string_path& path,
                                  std::shared_ptr<TRandom> r);

  virtual photon generate(const beam&, const target&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return E_range_.width(); }

protected:
  double intensity(const double E, const double beam, const double vz) const;

private:
  const realistic_target target_;  // target RL info
  const double E_beam_;            // (maximum) electron beam energy
  const interval<double> E_range_; // photon energy range
  const double max_;               // the maximum intensity.
};

// virtual photons
class vphoton : public photon_generator {
public:
  vphoton(const configuration& cf, const string_path& path,
          std::shared_ptr<TRandom> r);

  virtual photon generate(const beam&, const target&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const {
    return logy_range_.width() * logQ2_range_.width();
  }

protected:
  double flux(const double Q2, const double y, const particle& beam,
              const particle& target) const {
    return physics::gamma_t_log(Q2, y, beam, target);
  }

private:
  double calc_max_flux(const configuration& cf) const;
  interval<double> calc_max_Q2_range(const configuration& cf) const;
  interval<double> calc_max_W2_range(const configuration& cf) const;

  // primary kinematic boundaries
  const interval<double> y_range_;
  const interval<double> Q2_range_;
  // derived kinematic boundaries
  const interval<double> logy_range_;
  const interval<double> logQ2_range_;
  // additional cuts
  const interval<double> W2_range_;

  // maximum flux
  const double max_;
};

} // namespace initial
} // namespace lager

#endif
