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

#include "photon_gen.hh"
#include <TF1.h>
#include <TF2.h>
#include <TSpline.h>
#include <cmath>
#include <lager/core/assert.hh>
#include <lager/physics/kinematics.hh>
#include <lager/physics/photon.hh>

// =======================================================================================
// local utility functions
// =======================================================================================
namespace {
const lager::translation_map<lager::initial::bremsstrahlung::model>
    bs_model_translator{
        {"flat", lager::initial::bremsstrahlung::model::FLAT},
        {"param", lager::initial::bremsstrahlung::model::PARAM},
        {"approx", lager::initial::bremsstrahlung::model::APPROX},
        {"exact", lager::initial::bremsstrahlung::model::EXACT}};
// Bremsstrahlung spectrum for a 1%, 5% and 10% r.l. radiator interpolated
// between the  resultes of the exact calculation by Tsai and Whitis
// (SLAC-PUB-184 1966  (Table I))
inline double bremsstrahlung_intensity_001_param(const double E0,
                                                 const double k) {
  static double xv[] = {0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.82,
                        0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1};
  static double yv[] = {1.20417, 1.10071, 1.01740, 0.95406, 0.91055, 0.88668,
                        0.88222, 0.89671, 0.90181, 0.90760, 0.91407, 0.92119,
                        0.92889, 0.93710, 0.94566, 0.95421, 0.96164, 0.95451};
  static const TSpline3 brems10{"brems001", xv, yv, 18};
  return brems10.Eval(k / E0) * 0.01 / k;
}
inline double bremsstrahlung_intensity_005_param(const double E0,
                                                 const double k) {
  static double xv[] = {0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.82,
                        0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1};
  static double yv[] = {1.18110, 1.07713, 0.99387, 0.93046, 0.88616, 0.86013,
                        0.85123, 0.85735, 0.86001, 0.86302, 0.86625, 0.86956,
                        0.87272, 0.87534, 0.87674, 0.87538, 0.86653, 0.79765};
  static const TSpline3 brems10{"brems005", xv, yv, 18};
  return brems10.Eval(k / E0) * 0.05 / k;
}
inline double bremsstrahlung_intensity_010_param(const double E0,
                                                 const double k) {
  static double xv[] = {0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.82,
                        0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1};
  static double yv[] = {1.15277, 1.04815, 0.96490, 0.90130, 0.85591, 0.82716,
                        0.81289, 0.80929, 0.80921, 0.80908, 0.80873, 0.80789,
                        0.80616, 0.80289, 0.79691, 0.78565, 0.76136, 0.64338};
  static const TSpline3 brems10{"brems010", xv, yv, 18};
  return brems10.Eval(k / E0) * 0.1 / k;
}
// exact BS spectrum from Tsai and Whitis, gives the values consistent with the
// interpolation results above
inline double bremsstrahlung_intensity_exact(const double rad_len,
                                             const double E0, const double k) {
  static TF1 integrand(
      "I_g_gen1_integrand",
      [](double* ttprime, double* uu) {
        const double tprime = ttprime[0];
        const double u = uu[0]; // u = k/E0
        const double ln_u_inv = std::log(1 / u);
        const double a = 4. * tprime / 3.;
        const double fact0 = std::exp(7. * tprime / 9.) / tgamma(a + 1);
        const double fact1 = std::pow(ln_u_inv, a);
        const double term0 = u;
        double term1 = 0;
        // first 10 terms of perturbative expansion in eq 24
        for (int i = 0; i < 10; ++i) {
          double val = 1. / (tgamma(i + 1) * (i + a + 1));
          val *= ((4. / 3.) * std::pow(-1, i)) - u * u;
          val *= std::pow(ln_u_inv, i + 1);
          term1 += val;
        }
        return fact0 * fact1 * (term0 + term1);
      },
      0, 2, 1);
  integrand.SetParameter(0, k / E0);
  const double integral = integrand.Integral(0, rad_len);
  return std::exp(-7. * rad_len / 9.) * integral / k;
}
} // namespace

namespace lager {
namespace initial {

// =======================================================================================
// no-photon implementation
// =======================================================================================
no_photon::no_photon(const configuration& cf, const string_path& path,
                     std::shared_ptr<TRandom> r)
    : photon_generator{std::move(r)} {}
photon no_photon::generate(const beam&, const target&) { return {1.}; }

// =======================================================================================
// bremsstrahlung constructor
// =======================================================================================
bremsstrahlung::bremsstrahlung(const configuration& cf, const string_path& path,
                               std::shared_ptr<TRandom> r)
    : photon_generator{std::move(r)}
    , model_{cf.get<bremsstrahlung::model>(path / "model", bs_model_translator)}
    , rl_{(model_ != model::FLAT) ? cf.get<double>(path / "rl") : -1}
    , E_beam_{cf.get<double>("beam/lepton/energy")}
    , E_range_{cf.get_range<double>(path / "E_range")}
    , max_{intensity(E_range_.min, E_beam_)} {
  // initial info
  LOG_INFO("bremsstrahlung", "Maximum primary electron beam energy [GeV]: " +
                                 std::to_string(E_beam_));
  LOG_INFO("bremsstrahlung",
           "Photon energy range [GeV]: " + std::to_string(E_range_.min) + ", " +
               std::to_string(E_range_.max) + "]");
  // validate the setup
  tassert(E_range_.max <= E_beam_,
          "Photon energy cannot exceed electron beam energy");
  tassert(E_range_.width() > 0,
          "Ensure Emin < Emax for the photon energy range");
  tassert(E_range_.min > 0, "Ensure Emin > 0 for the photon beam energy");
  // Radiation length and model info
  if (model_ == model::FLAT) {
    LOG_INFO("bremsstrahlung", "Using a flat BS distribution.");
  } else if (model_ == model::PARAM) {
    LOG_INFO("bremsstrahlung",
             "Using parameterizatoin of the exact BS spectrum.");
    LOG_INFO("bremsstrahlung", "RL: " + std::to_string(rl_));
    if (rl_ != .01 && rl_ != .05 && rl_ != .10) {
      LOG_ERROR(
          "bremsstrahlung",
          "No parameterization available for the requested radiation length.");
      LOG_ERROR("bremsstrahlung", "Available RL values: [.01, .05, .10].");
      LOG_ERROR("bremsstrahlung", "Use 'approx' instead of 'param' to use an "
                                  "arbitrary radiation length.");
      throw cf.value_error("rl");
    }
  } else if (model_ == model::EXACT) {
    LOG_INFO("bremsstrahlung", "Using exact of BS spectrum.");
    LOG_INFO("bremsstrahlung", "RL: " + std::to_string(rl_));
  } else { // APPROX
    LOG_INFO("bremsstrahlung", "Using approximation of the exact BS spectrum");
    LOG_INFO("bremsstrahlung", "RL: " + std::to_string(rl_));
  }
}

// =======================================================================================
// Generate a new bremstrahlung photon
// =======================================================================================
photon bremsstrahlung::generate(const beam& lepton, const target& targ) {
  tassert(lepton.particle().energy() > E_range_.min,
          "Beam energy below bremsstrahlung generation window");
  tassert(lepton.particle().energy() <= E_beam_,
          "Beam energy higher than maximum electron beam energy.");
  // generate a value for E
  const double E = rng()->Uniform(E_range_.min, E_range_.max);
  LOG_JUNK("bremsstrahlung", "Generated E: " + std::to_string(E));

  // check if this value is in the allowed range
  const interval<double> Elim{E_range_.min,
                              fmin(lepton.particle().energy(), E_range_.max)};
  if (Elim.excludes(E)) {
    LOG_JUNK("bremsstrahlung", "Value outside the valid E range");
    return photon{0.};
  }

  photon pd = photon::make_real(lepton.particle(), targ.particle(), E,
                                intensity(E, lepton.particle().energy()));

  // that's all!
  return pd;
}

// =======================================================================================
// return the bremstrahlung intensity for the prefered parameterization
// =======================================================================================
double bremsstrahlung::intensity(const double E, const double E_beam) const {
  if (model_ == model::FLAT) {
    return 1.;
  } else if (model_ == model::PARAM) {
    if (rl_ == .01) {
      return bremsstrahlung_intensity_001_param(E_beam, E);
    } else if (rl_ == .05) {
      return bremsstrahlung_intensity_005_param(E_beam, E);
    } else if (rl_ == .10) {
      return bremsstrahlung_intensity_010_param(E_beam, E);
    } else {
      return 0.; // can never happen
    }
  } else if (model_ == model::EXACT) {
    return bremsstrahlung_intensity_exact(rl_, E_beam, E);
  } else { // APPROX
    return physics::bremsstrahlung_approx(rl_, E_beam, E);
  }
}
// =======================================================================================
// bremsstrahlung_realistic_target constructor
// =======================================================================================
bremsstrahlung_realistic_target::bremsstrahlung_realistic_target(
    const configuration& cf, const string_path& path,
    std::shared_ptr<TRandom> r)
    : photon_generator{std::move(r)}
    , target_{cf, path}
    , E_beam_{cf.get<double>("beam/lepton/energy")}
    , E_range_{cf.get_range<double>(path / "E_range")}
    , max_{intensity(E_range_.min, E_beam_, target_.back())} {
  // initial info
  LOG_INFO("bremsstrahlung_realistic_target",
           "Maximum primary electron beam energy [GeV]: " +
               std::to_string(E_beam_));
  LOG_INFO("bremsstrahlung_realistic_target",
           "Photon energy range [GeV]: " + std::to_string(E_range_.min) + ", " +
               std::to_string(E_range_.max) + "]");
  // validate the setup
  tassert(E_range_.max <= E_beam_,
          "Photon energy cannot exceed electron beam energy");
  tassert(E_range_.width() > 0,
          "Ensure Emin < Emax for the photon energy range");
  tassert(E_range_.min > 0, "Ensure Emin > 0 for the photon beam energy");
  tassert(target_.length() >= 0,
          "Ensure min <= max for the target z-coordinate range");
}

// =======================================================================================
// Generate a new bremstrahlung photon
// =======================================================================================
photon bremsstrahlung_realistic_target::generate(const beam& lepton,
                                                 const target& targ) {
  tassert(lepton.particle().energy() > E_range_.min,
          "Beam energy below bremsstrahlung generation window");
  tassert(lepton.particle().energy() <= E_beam_,
          "Beam energy higher than maximum electron beam energy.");
  // generate a value for E
  const double E = rng()->Uniform(E_range_.min, E_range_.max);
  LOG_JUNK("bremsstrahlung_realistic_target",
           "Generated E: " + std::to_string(E));

  // check if this value is in the allowed range
  const interval<double> Elim{E_range_.min,
                              fmin(lepton.particle().energy(), E_range_.max)};
  if (Elim.excludes(E)) {
    LOG_JUNK("bremsstrahlung", "Value outside the valid E range");
    return photon{0.};
  }

  photon pd = photon::make_real(
      lepton.particle(), targ.particle(), E,
      intensity(E, lepton.particle().energy(), lepton.particle().vertex().z()));

  // that's all!
  return pd;
}

// =======================================================================================
// return the bremstrahlung intensity for the prefered parameterization
// =======================================================================================
double bremsstrahlung_realistic_target::intensity(const double E,
                                                  const double E_beam,
                                                  const double vz) const {
  return bremsstrahlung_intensity_exact(target_.total_rl(vz), E_beam, E);
}

// =======================================================================================
// constructor for vphoton
// =======================================================================================
vphoton::vphoton(const configuration& cf, const string_path& path,
                 std::shared_ptr<TRandom> r)
    : photon_generator{std::move(r)}
    , y_range_{cf.get_range<double>(path / "y_range")}
    , logy_range_{std::log(y_range_.min), std::log(y_range_.max)}
    , Q2_range_{calc_max_Q2_range(cf)}
    , logQ2_range_{std::log(Q2_range_.min), std::log(Q2_range_.max)}
    , W2_range_{calc_max_W2_range(cf)}
    , max_{calc_max_flux(cf)} {
  // initial info
  LOG_INFO("vphoton", "Q2 range [GeV^2]: [" + std::to_string(Q2_range_.min) +
                          ", " + std::to_string(Q2_range_.max) + "]");
  LOG_INFO("vphoton", "y range [GeV^2]: [" + std::to_string(y_range_.min) +
                          ", " + std::to_string(y_range_.max) + "]");
  if (cf.get_optional_range<double>("W_range")) {
    LOG_INFO("vphoton", "W range (optional) [GeV}: [" +
                            std::to_string(sqrt(W2_range_.min)) + ", " +
                            std::to_string(sqrt(W2_range_.max)) + "]");
  } else {
    LOG_INFO("vphoton", "W range (default) [GeV}: [" +
                            std::to_string(sqrt(W2_range_.min)) + ", " +
                            std::to_string(sqrt(W2_range_.max)) + "]");
  }
  // validate the setup
  tassert(y_range_.min > 0, "Ensure ymin > 0");
  tassert(y_range_.max <= 1, "Ensure ymax <= 1");
  tassert(Q2_range_.width() > 0,
          "Ensure Q2min < Q2max for the virtual photon Q2 range");
  tassert(y_range_.width() > 0,
          "Ensure ymin < ymax for the virtual photon y range");
}

// =======================================================================================
// generate a new virtual photon
// =======================================================================================
photon vphoton::generate(const beam& lepton, const target& targ) {

  // generate a value for Q2 and y
  const double y = exp(rng()->Uniform(logy_range_.min, logy_range_.max));
  const double Q2 = exp(rng()->Uniform(logQ2_range_.min, logQ2_range_.max));

  LOG_JUNK("vphoton",
           "Generated y: " + std::to_string(y) + " Q2: " + std::to_string(Q2));

  // check Q2 boundaries for validity
  if (physics::Q2_range(lepton.particle(), targ.particle(), y).excludes(Q2)) {
    LOG_JUNK("vphoton", "Values outside of valid Q2 range");
    return {0.};
  }

  photon pd =
      photon::make_virtual(lepton.particle(), targ.particle(), Q2, y,
                           flux(Q2, y, lepton.particle(), targ.particle()),
                           rng()->Uniform(0, TMath::TwoPi()));

  LOG_JUNK("vphoton", "nu: " + std::to_string(pd.nu()) +
                          " W2: " + std::to_string(pd.W2()) +
                          " x: " + std::to_string(pd.x()));

  // check if the invariants are in the range we want
  if (W2_range_.excludes(pd.W2())) {
    LOG_JUNK("vphoton", "Values outside of valid W2 range");
    pd.update_cross_section(0);
    return pd;
  }

  LOG_JUNK("vphoton", "xsec: " + std::to_string(pd.cross_section()) + "( < " +
                          std::to_string(max_) + ")");

  // bail if our cross section is not positive
  if (pd.cross_section() < 0) {
    pd.update_cross_section(0);
  }

  // that's all
  return pd;
}

// =======================================================================================
// calculate an upper limit for the flux for the requested kinematic limits
// =======================================================================================
double vphoton::calc_max_flux(const configuration& cf) const {
  const particle beam{cf.get<std::string>("beam/lepton/particle_type"),
                      cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                      cf.get<double>("beam/lepton/energy")};
  const particle target{cf.get<std::string>("beam/ion/particle_type"),
                        cf.get_vector3<particle::XYZVector>("beam/ion/dir"),
                        cf.get<double>("beam/ion/energy")};
  TF2 fflux(
      "flux",
      [=](double* Q2y, double* par = 0x0) {
        return this->flux(Q2y[0], Q2y[1], beam, target);
      },
      Q2_range_.min, Q2_range_.max, y_range_.min, y_range_.max, 0);
  double Q2, y;
  fflux.GetMaximumXY(Q2, y);
  return flux(Q2, y, beam, target) * 1.01;
}

// =======================================================================================
// calculate max Q2 range, assuming the y-range is already set
// =======================================================================================
interval<double> vphoton::calc_max_Q2_range(const configuration& cf) const {
  // get beam and target info
  const particle beam{cf.get<std::string>("beam/lepton/particle_type"),
                      cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                      cf.get<double>("beam/lepton/energy")};
  const particle target{cf.get<std::string>("beam/ion/particle_type"),
                        cf.get_vector3<particle::XYZVector>("beam/ion/dir"),
                        cf.get<double>("beam/ion/energy")};

  // Cap the minimum Q2 at lepton mass squared
  const double Q2min =
      fmax(physics::Q2_range(beam, target, y_range_.min).min, 1.0e-12);

  // In general, the maximum Q2 range is determined through the kinematics
  // of t-channel scattering (falling function of y), as well as the
  // requirement that the final state contain at least te target mass (W2min
  // = M2_target; a rising function of y).
  //
  // Q2 max is reached for the point where both values cross,  which is (up
  // to a precision of lepton mass squared): y = 2E/(M+2E).
  //
  // alternatively, if a y-cut is set below this point, the maximum is
  // reached at maximum y
  //
  // finally, we also allow for the user to manually set a Q2 range
  const double E = (beam.p()).Dot(target.p()) / target.mass();
  // get Q2 range from kinematics
  const double Q2max =
      physics::Q2_range(beam, target,
                        fmin(2. * E / (target.mass() + 2. * E), y_range_.max))
          .max;
  // check if there is an alternative user-defined range
  const auto opt_range_ = cf.get_optional_range<double>("photon/Q2_range");
  if (opt_range_) {
    // only return a kinematically allowed range
    return {fmax(opt_range_->min, Q2min), fmin(opt_range_->max, Q2max)};
  }
  return {Q2min, Q2max};
}

// =======================================================================================
// allow for the user to set a W2 range
// =======================================================================================
interval<double> vphoton::calc_max_W2_range(const configuration& cf) const {
  // get beam and target info
  const particle beam{cf.get<std::string>("beam/lepton/particle_type"),
                      cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                      cf.get<double>("beam/lepton/energy")};
  const particle target{cf.get<std::string>("beam/ion/particle_type"),
                        cf.get_vector3<particle::XYZVector>("beam/ion/dir"),
                        cf.get<double>("beam/ion/energy")};

  // at least target in final state, at most s in final state
  const double W2min = target.mass2();
  const double W2max = (beam.p() + target.p()).M2();

  // check if the user requested a different range range
  const auto opt_range_ = cf.get_optional_range<double>("photon/W_range");

  if (opt_range_) {
    // ensure the user-set range is kinematically allowed
    return {fmax(opt_range_->min * opt_range_->min, W2min),
            fmin(opt_range_->max * opt_range_->max, W2max)};
  }
  return {W2min, W2max};
}

} // namespace initial
} // namespace lager
