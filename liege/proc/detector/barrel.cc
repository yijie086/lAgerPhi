#include "barrel.hh"

namespace liege {
namespace detector {

barrel::barrel(const configuration& cf, const string_path& path,
               std::shared_ptr<TRandom> r)
    : barrel::base_type{r}
    , name_{cf.get<std::string>(path / "name")}
    , id_{cf.get<int>(path / "id", 0)}
    , theta_{cf.get_range<double>(path / "acceptance" / "theta") *
             TMath::DegToRad()}
    , p_{cf.get_range<double>(path / "acceptance" / "p")}
    , acceptance_{cf.get<double>(path / "acceptance" / "acceptance")}
    , pid_{cf.get_vector<int>(path / "acceptance" / "pid")}
    , p_smear_{cf.get<double>(path / "smearing" / "p", 0.)}
    , theta_smear_{cf.get<double>(path / "smearing" / "theta", 0.) / 1000.}
    , phi_smear_{cf.get<double>(path / "smearing" / "phi", 0.) / 1000.} {
  LOG_INFO(name_, "ID: " + std::to_string(id_));
  LOG_INFO(name_, "Theta range [rad]: [" + std::to_string(theta_.min) + ", " +
                      std::to_string(theta_.max) + "]");
  LOG_INFO(name_, "p range [GeV]: [" + std::to_string(p_.min) + ", " +
                      std::to_string(p_.max) + "]");
  LOG_INFO(name_, "Acceptance: " + std::to_string(acceptance_));
  LOG_INFO(name_, "PID: " + stringify(pid_));
  if (p_smear_ > 0 || theta_smear_ > 0 || phi_smear_ > 0) {
    LOG_INFO(name_, "Momentum smearing: " + std::to_string(p_smear_));
    LOG_INFO(name_, "Theta smearing [rad]: " + std::to_string(theta_smear_));
    LOG_INFO(name_, "Phi smearing [rad]: " + std::to_string(phi_smear_));
  } else {
    LOG_INFO(name_, "No smearing");
  }
}

ROOT::Math::PxPyPzMVector barrel::detected_track(const particle& part) const {
  const double p =
      (p_smear_ > 0) ? rng()->Gaus(part.momentum(), p_smear_ * part.momentum())
                     : part.momentum();
  const double theta = (theta_smear_ > 0)
                           ? rng()->Gaus(part.theta(), theta_smear_)
                           : part.theta();
  const double phi =
      (phi_smear_ > 0) ? rng()->Gaus(part.phi(), phi_smear_) : part.phi();
  const double px = p * sin(theta) * cos(phi);
  const double py = p * sin(theta) * sin(phi);
  const double pz = p * cos(theta);
  TVector3 det_vec{px, py, pz};
  if (p_smear_ > 0 && theta_smear_ > 0 && phi_smear_ > 0) {
    LOG_JUNK2(name_, "Smeared variables (P, theta, phi): (" +
                         std::to_string(part.momentum()) + ", " +
                         std::to_string(part.theta()) + ", " +
                         std::to_string(part.phi()) + ") --> (" +
                         std::to_string(p) + ", " + std::to_string(theta) +
                         ", " + std::to_string(phi) + ")")
  }
  return {det_vec.X(), det_vec.Y(), det_vec.Z(), part.mass()};
}

void barrel::process(event& e) const {
  for (auto& part : e) {
    if (part.final_state() &&
        std::any_of(pid_.begin(), pid_.end(), [&](const auto& good_pid) {
          return part.type<int>() == good_pid;
        })) {
      LOG_JUNK2(name_, "Found matching final state particle " + part.name() +
                           " (status: " + std::to_string(part.status<int>()) +
                           ", momentum: " + std::to_string(part.momentum()) +
                           ")");
      // check barrel cuts
      if (p_.includes(part.momentum())) {
        LOG_JUNK2(name_, "Momentum cut for " + part.name() + ": SUCCESS");
        LOG_JUNK2(name_, "True theta: " + std::to_string(part.theta()) +
                             ", phi: " + std::to_string(part.phi()));
        if (theta_.includes(part.theta())) {
          LOG_JUNK2(name_, "Barrel cut for " + part.name() + ": SUCCESS");
          if (acceptance_ == 1. || rng()->Uniform(0, 1.) < acceptance_) {
            LOG_JUNK2(name_,
                      "Flat acceptance for " + part.name() + ": SUCCESS");
            auto detected = detected_track(part);
            e.add_detected(
                {part,
                 {detected.X(), detected.Y(), detected.Z(), detected.E()},
                 id_});
          } else {
            LOG_JUNK2(name_, "Flat acceptance for " + part.name() + ": FAILED");
          }
        } else {
          LOG_JUNK2(name_, "Barrel cut for " + part.name() + ": FAILED");
        }
      } else {
        LOG_JUNK2(name_, "Momentum cut for " + part.name() + ": FAILED");
      }
    }
  }
}
} // namespace detector
} // namespace liege
