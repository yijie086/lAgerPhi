#include "bremsstrahlung.hh"
#include <TLorentzVector.h>
#include <pcsim/core/assert.hh>
#include <pcsim/core/logger.hh>

namespace pcsim {
namespace physics {

// bremsstrahlung constructor
bremsstrahlung::bremsstrahlung(const ptree& settings, const string_path& path)
    : configurable{settings, path}
    , model_{get_model()}
    , rl_{(model_ == model::APPROX) ? conf().get<double>("rl") : 0.}
    , E0_{conf().get<double>("electron_energy")}
    , range_{conf().get_range<double>("range")}
    , max_{calc_max()} {
  tassert(range_.max <= E0_, "Photon energy cannot exceed beam energy");
  tassert(range_.width() > 0, "Ensure Emin < Emax for the photon beam energy");
  tassert(range_.min > 0, "Ensure Emin > 0 for the photon beam energy");
  LOG_INFO("bremsstrahlung",
           "Primary electron beam energy: " + std::to_string(E0_) + " GeV");
  LOG_INFO("bremsstrahlung", "Photon energy range: [" +
                                 std::to_string(range_.min) + ", " +
                                 std::to_string(range_.max) + "] GeV");
}

interval<double> bremsstrahlung::calc_s_range() const {
  interval<double> E_range = conf().get_range<double>("range");
  const TLorentzVector photon_min{0, 0, E_range.min, E_range.min};
  const TLorentzVector photon_max{0, 0, E_range.max, E_range.max};
  const TLorentzVector target{0, 0, 0, physics::M_PROTON};
  return {(target + photon_min).M2(), (target + photon_max).M2()};
}

// private utility function to read the correct BS model from the configuration
bremsstrahlung::model bremsstrahlung::get_model() const {
  if (conf().get<std::string>("type") == "flat") {
    LOG_INFO("bremsstrahlung", "Using a flat BS distribution");
    return model::FLAT;
  } else if (conf().get<std::string>("type") == "param") {
    LOG_INFO("bremsstrahlung",
             "Using parameterization of the exact BS spectrum for 0.1 RL.");
    return model::PARAM;
  } else {
    LOG_INFO("bremsstrahlung", "RL: " + conf().get<std::string>("rl"));
    LOG_INFO("bremsstrahlung", "Using approximation of the exact BS spectrum)");
    return model::APPROX;
  }
}

// private utility function to calculate the maximum intensity needed for
// accept-reject, based on the chosen model
double bremsstrahlung::calc_max() const {
  double max = 0;
  const double jacobian = physics::ONE_OVER_2M_PROTON;
  // for a flat distribution, the max is given by 1/range
  if (model_ == model::FLAT) {
    max = 1. / jacobian;
    // for the parameterization, the spectrum is monotonously falling
  } else if (model_ == model::PARAM) {
    max = physics::bremsstrahlung_intensity_10_param(E0_, range_.min);
    // the approximate spectrum is also monotonously falling
  } else {
    max = physics::bremsstrahlung_intensity_approx(rl_, E0_, range_.min);
  }
  return max * jacobian;
}

} // physics
} // pcsim
