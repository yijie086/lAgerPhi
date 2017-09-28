#include "pc.hh"
#include <pcsim/core/assert.hh>
#include <pcsim/core/logger.hh>
#include <pcsim/physics/constants.hh>
#include <pcsim/core/accept_reject.hh>

namespace pcsim {
namespace physics {

// =============================================================================
// Pc decay function object implementation
// =============================================================================

// constructor
pc_decay::pc_decay(const configuration& cf, const string_path& path,
                   std::shared_ptr<TRandom> rng)
    : configurable{cf, path}, mode_{get_mode()}, rng_{std::move(rng)} {}
// constructor utility function
pc_decay::mode pc_decay::get_mode() const {
  if(conf().get<std::string>("type") == "iso") {
    LOG_INFO("pc_decay", "Using isotropical Pc decay distribution");
    return mode::ISO;
  } else if (conf().get<std::string>("type") == "5/2+") {
    LOG_INFO("pc_decay", "Using fit to 5/2+ decay distribution");
    return mode::S52_PLUS;
  } else if (conf().get<std::string>("type") == "5/2-") {
    LOG_INFO("pc_decay", "Using fit to 5/2- decay distribution");
    return mode::S52_MINUS;
  } else if (conf().get<std::string>("type") == "3/2+") {
    LOG_INFO("pc_decay", "Using fit to 3/2+ decay distribution");
    return mode::S32_PLUS;
  } else if (conf().get<std::string>("type") == "3/2-") {
    LOG_INFO("pc_decay", "Using fit to 3/2- decay distribution");
    return mode::S32_MINUS;
  } else {
    LOG_ERROR("pc_decay", "Invalid spin/parity combination");
    throw conf().value_error("type");
  }
}

// Simulate the Pc --> J/Psi,p decay, using the appropriate angular
// distribution for the requested spin/parity mode
void pc_decay::operator()(const TLorentzVector& pc, TLorentzVector& proton,
                          TLorentzVector& jpsi) const {
  // work in the Pc CM frame
  const double phi = rng_->Uniform(0., TMath::TwoPi());
  double ctheta = -1;
  if (mode_ == mode::ISO) {
    ctheta = rng_->Uniform(-1, 1);
  } else if (mode_ == mode::S52_PLUS) {
    // result from a pol6 fit to a digitized version of figure 6c from
    // PRD92-034022(2015)
    ctheta = accept_reject_1D(rng_, {-1, 1},
                              [](const double x) {
                                const double x2 = x * x;
                                const double x3 = x2 * x;
                                const double x4 = x3 * x;
                                const double x5 = x4 * x;
                                const double x6 = x5 * x;
                                return .149211 - 0.194418 * x - 0.563191 * x2 +
                                       0.374024 * x3 + 0.658942 * x4 +
                                       0.110057 * x5 + 0.0931712 * x6;
                              },
                              0.63);
  } else if (mode_ == mode::S52_MINUS) {
    // result from a pol7 fit to a digitized version of figure 5c from
    // PRD92-034022(2015)
    ctheta = accept_reject_1D(rng_, {-1, 1},
                              [](const double x) {
                                const double x2 = x * x;
                                const double x3 = x2 * x;
                                const double x4 = x3 * x;
                                const double x5 = x4 * x;
                                const double x6 = x5 * x;
                                const double x7 = x6 * x;
                                return 1.31241 - 1.19802 * x + 1.58351 * x2 +
                                       17.1514 * x3 + 20.8306 * x4 -
                                       4.43848 * x5 + 2.67151 * x6 +
                                       6.06378 * x7;
                              },
                              44.06);
  } else if (mode_ == mode::S32_PLUS) {
    // result from a expo fit to a digitized version of figure 5b from
    // PRD92-034022(2015)
    ctheta = accept_reject_1D(
        rng_, {-1, 1}, [](const double x) { return exp(-5.944 - x); }, 0.00713);
  } else if (mode_ == mode::S32_MINUS) {
    // result from a pol2 fit to a digitized version of figure 6b from
    // PRD92-034022(2015)
    ctheta =
        accept_reject_1D(rng_, {-1, 1},
                         [](const double x) {
                           const double x2 = x * x;
                           return 0.00845846 - 0.0128146 * x + 0.00526053 * x2;
                         },
                         0.0266);
  } else {
    tassert(false, "SHOULD NEVER HAPPEN");
  }
  const double theta = acos(ctheta);
  const double E_cm = pc.M();
  const double Ep = (E_cm * E_cm + M2_PROTON - M2_JPSI) / (2 * E_cm);
  const double Ej = (E_cm * E_cm - M2_PROTON + M2_JPSI) / (2 * E_cm);
  const double Pp = sqrt(Ep * Ep - M2_PROTON);
  const double Pj = sqrt(Ej * Ej - M2_JPSI);
  TVector3 mom;
  // J/Psi going forward, recoil going backward
  mom.SetMagThetaPhi(Pj, theta, phi);
  jpsi = {mom, Ej};
  mom.SetMagThetaPhi(Pp, theta, phi);
  proton = {-mom, Ep};
  // boost to lab frame
  const auto beta = pc.BoostVector();
  proton.Boost(beta);
  jpsi.Boost(beta);
  // all done!
}

} // physics
} // pcsim
