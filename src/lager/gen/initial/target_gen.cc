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

#include "target_gen.hh"

#include <lager/physics/fermi.hh>

#include <cmath>

#include <TF1.h>

namespace {
// calculate A based on the ion mass
int calc_A(const std::string& ion_name) {
  lager::particle ion{ion_name};
  return std::nearbyint(ion.mass() / .938272);
}
} // namespace

namespace lager {
namespace initial {

// estimated target info from config file to be used by lA generators to
// estimate phase-space and cross section ranges
particle estimated_target(const configuration& cf) {
  const particle ion{cf.get<std::string>("beam/ion/particle_type"),
                     cf.get_vector3<particle::XYZVector>("beam/ion/dir"),
                     cf.get<double>("beam/ion/energy")};
  if (cf.get<std::string>("target/type") == "primary") {
    return ion;
    // right now only other option is fermi momentum. In that case scale target
    // mass to proton mass and add extra 1.2GeV in fermi-momentum
  } else {
    particle proton{pdg_id::p,
                    cf.get_vector3<particle::XYZVector>("beam/ion/dir"),
                    sqrt(.938372 * .938272 + 1.2 * 1.2)};
    lager::particle::Boost boost_to_lab{-ion.p().BoostToCM()};
    proton.boost(boost_to_lab);
    return proton;
  }
}

// =======================================================================================
// target beam constructor
// =======================================================================================
primary_target::primary_target(const configuration& cf, const string_path& path,
                               std::shared_ptr<TRandom> r)
    : target_generator{std::move(r)}
    , target_{cf.get<std::string>("beam/ion/particle_type"),
              particle::status_code::SECONDARY_BEAM} {
  LOG_INFO("primary_target", "Using primary ion beam as target");
}
target primary_target::generate(const beam& ion) {
  return target::make_target(ion.particle(), target_);
}

// =======================================================================================
// fermi87 target
// =======================================================================================
fermi87::fermi87(const configuration& cf, const string_path& path,
                 std::shared_ptr<TRandom> r)
    : target_generator{std::move(r)}
    , A_{calc_A(cf.get<std::string>("beam/ion/particle_type"))}
    , nucleon_{cf.get<std::string>(path / "nucleon")}
    , k_max_{cf.get<double>(path / "k_max")}
    , norm_{calc_norm()}
    , xs_max_{calc_max()} {
  LOG_INFO("initial::fermi87",
           "Nucleus: " + cf.get<std::string>("beam/ion/particle_type"));
  LOG_INFO("initial::fermi87", "Calculated A: " + std::to_string(A_));
  LOG_INFO("initial::fermi87", "Nucleon: " + nucleon_.name());
  LOG_INFO("initial::fermi87", "k_max [GeV]: " + std::to_string(k_max_));
  LOG_DEBUG("initial::fermi87",
            "Normalization factor: " + std::to_string(norm_));
  LOG_DEBUG("initial::fermi87", "P.D.F. max:" + std::to_string(xs_max_));
}
target fermi87::generate(const beam& ion) {
  // Generate a nucleon momentum, theta and phi
  const double P = rand_f(
      {0., k_max_}, [&](const double P) { return pdf(P); }, xs_max_ * 1.01);
  const double theta = acos(rng()->Uniform(-1, 1));
  const double phi = rng()->Uniform(0., TMath::TwoPi());
  LOG_DEBUG("initial::fermi87",
            "Selected a target nucleon with (P, theta, phi): (" +
                std::to_string(P) + ", " +
                std::to_string(theta * TMath::RadToDeg()) + ", " +
                std::to_string(phi * TMath::RadToDeg()) + ")");
  const particle::Polar3DVector p3{P, theta, phi};
  return target::make_target(ion.particle(),
                             {nucleon_.type(),
                              {p3.X(), p3.Y(), p3.Z()},
                              particle::status_code::SECONDARY_BEAM});
}

// fermi87 private
double fermi87::calc_norm() const {
  TF1 f{"tmp", [=](double* pp, double*) { return physics::fermi87(pp[0], A_); },
        0, k_max_, 0};
  double integral = f.Integral(0, k_max_);
  tassert(integral > 0,
          "Fermi momentum normalization should be positive number");
  return 1. / integral;
}
double fermi87::pdf(const double P) const {
  return norm_ * physics::fermi87(P, A_);
}
// maximum is reached for nucleons at rest
double fermi87::calc_max() const {
  TF1 f{"tmp", [=](double* pp, double*) { return pdf(pp[0]); }, 0, k_max_, 0};
  return f.GetMaximum(0, k_max_);
}

} // namespace initial
} // namespace lager
