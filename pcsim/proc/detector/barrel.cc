#include "barrel.hh"

namespace pcsim {
namespace detector {

barrel::barrel(const configuration& cf, const string_path& path,
                           std::shared_ptr<TRandom> r)
    : barrel::base_type{r}
    , name_{cf.get<std::string>(path / "name")}
    , id_{cf.get<int>(path / "id", 0)}
    , theta_{cf.get_range<double>(path / "acceptance" / "theta") * TMath::DegToRad()}
    , p_{cf.get_range<double>(path / "acceptance" / "p")}
    , acceptance_{cf.get<double>(path / "acceptance" / "acceptance")}
    , pid_{cf.get_vector<int>(path / "acceptance" / "pid")}
    , p_smear_{cf.get<double>(path / "smearing" / "p", 0.)}
    , theta_smear_{cf.get<double>(path / "smearing" / "theta", 0.) / 1000.} {
  LOG_INFO(name_, "ID: " + std::to_string(id_));
  LOG_INFO(name_, "Theta range [rad]: [" +
                      std::to_string(theta_.min) + ", " +
                      std::to_string(theta_.max) + "]");
  LOG_INFO(name_, "p range [GeV]: [" + std::to_string(p_.min) + ", " +
                      std::to_string(p_.max) + "]");
  LOG_INFO(name_, "Acceptance: " + std::to_string(acceptance_));
  LOG_INFO(name_, "PID: " + stringify(pid_));
  if (p_smear_ > 0 || theta_smear_ > 0) {
    LOG_INFO(name_, "Momentum smearing: " + std::to_string(p_smear_));
    LOG_INFO(name_,
             "Angle smearing [rad]: " + std::to_string(theta_smear_));
  } else {
    LOG_INFO(name_, "No smearing");
  }
}

ROOT::Math::PxPyPzMVector
barrel::detected_track(const particle& part, const double th_in,
                             const double th_out) const {
  const double p =
      (p_smear_ > 0) ? rng()->Gaus(part.momentum(), p_smear_ * part.momentum())
                     : part.momentum();
  const double thx =
      (th_in_smear_ > 0) ? rng()->Gaus(th_in, th_in_smear_) : th_in;
  const double thy =
      (th_out_smear_ > 0) ? rng()->Gaus(th_out, th_out_smear_) : th_out;
  const double px = p * sin(thx);
  const double py = p * sin(thy);
  const double pz = sqrt(p * p - px * px - py * py);
  TVector3 det_vec{px, py, pz};
  // undo rotation from track_th_in_out
  det_vec.RotateY(theta0_);
  det_vec.RotateZ(phi0_);
  if (p_smear_ > 0 && th_in_smear_ > 0 && th_out_smear_ > 0) {
    LOG_JUNK2(name_, "Smeared variables (P, th_in, th_out): (" +
                         std::to_string(part.momentum()) + ", " +
                         std::to_string(th_in) + ", " + std::to_string(th_out) +
                         ") --> (" + std::to_string(p) + ", " +
                         std::to_string(thx) + ", " + std::to_string(thy) + ")")
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
        auto [th_in, th_out, pz] = track_th_in_out_pz(part);
        LOG_JUNK2(name_, "True th_in: " + std::to_string(th_in) +
                             ", th_out: " + std::to_string(th_out));
        if (th_in_.includes(th_in) && th_out_.includes(th_out) && pz > 0) {
          LOG_JUNK2(name_, "Angular box cut for " + part.name() + ": SUCCESS");
          if (acceptance_ == 1. || rng()->Uniform(0, 1.) < acceptance_) {
            LOG_JUNK2(name_,
                      "Flat acceptance for " + part.name() + ": SUCCESS");
            auto detected = detected_track(part, th_in, th_out);
            e.add_detected(
                {part,
                 {detected.X(), detected.Y(), detected.Z(), detected.E()},
                 id_});
          } else {
            LOG_JUNK2(name_, "Flat acceptance for " + part.name() + ": FAILED");
          }
        } else {
          LOG_JUNK2(name_, "Angular box cut for " + part.name() + ": FAILED");
        }
      } else {
        LOG_JUNK2(name_, "Momentum cut for " + part.name() + ": FAILED");
      }
    }
  }
}
} // namespace detector
} // namespace pcsim
