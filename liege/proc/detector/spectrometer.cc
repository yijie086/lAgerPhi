#include "spectrometer.hh"

namespace liege {
namespace detector {

spectrometer::spectrometer(const configuration& cf, const string_path& path,
                           std::shared_ptr<TRandom> r)
    : spectrometer::base_type{r}
    , name_{cf.get<std::string>(path / "name")}
    , id_{cf.get<int>(path / "id", 0)}
    , theta0_{cf.get<double>(path / "position" / "theta0") * TMath::DegToRad()}
    , phi0_{cf.get<double>(path / "position" / "phi0") * TMath::DegToRad()}
    , p0_{cf.get<double>(path / "position" / "p0")}
    , acceptance_{cf.get<double>(path / "acceptance" / "acceptance")}
    , th_in_{cf.get_range<double>(path / "acceptance" / "th_in") / 1000.}
    , th_out_{cf.get_range<double>(path / "acceptance" / "th_out") / 1000.}
    , dp_{cf.get_range<double>(path / "acceptance" / "delta") * p0_ / 100.}
    , p_{dp_.min + p0_, dp_.max + p0_}
    , pid_{cf.get_vector<int>(path / "acceptance" / "pid")}
    , p_smear_{cf.get<double>(path / "smearing" / "p", 0.)}
    , th_in_smear_{cf.get<double>(path / "smearing" / "th_in", 0.) / 1000.}
    , th_out_smear_{cf.get<double>(path / "smearing" / "th_out", 0.) / 1000.} {
  LOG_INFO(name_, "ID: " + std::to_string(id_));
  LOG_INFO(name_,
           "Central angle [deg.]: " +
               std::to_string(cf.get<double>(path / "position" / "theta0")));
  LOG_INFO(name_, "Rotation [deg.]: " + std::to_string(cf.get<double>(
                                            path / "position" / "phi0")));
  LOG_INFO(name_, "Acceptance: " + std::to_string(acceptance_));
  LOG_INFO(name_, "Central momentum [GeV]: " + std::to_string(p0_));
  LOG_INFO(name_, "Inbending angle range [rad]: [" +
                      std::to_string(th_in_.min) + ", " +
                      std::to_string(th_in_.max) + "]");
  LOG_INFO(name_, "Outbending angle range [rad]: [" +
                      std::to_string(th_out_.min) + ", " +
                      std::to_string(th_out_.max) + "]");
  LOG_INFO(name_, "p range [GeV]: [" + std::to_string(p_.min) + ", " +
                      std::to_string(p_.max) + "]");
  LOG_INFO(name_, "PID: " + stringify(pid_));
  if (p_smear_ > 0 || th_in_smear_ > 0 || th_out_smear_ > 0) {
    LOG_INFO(name_, "Momentum smearing: " + std::to_string(p_smear_));
    LOG_INFO(name_,
             "Inbending angle smearing [rad]: " + std::to_string(th_in_smear_));
    LOG_INFO(name_, "Outbending angle smearing [rad]: " +
                        std::to_string(th_out_smear_));
  } else {
    LOG_INFO(name_, "No smearing");
  }
}

std::tuple<double, double, double>
spectrometer::track_th_in_out_pz(const particle& part) const {
  TVector3 track{part.p().X(), part.p().Y(), part.p().Z()};
  // rotate about Z-axis to coordinate system were X-axis is in the horizontal
  // plane
  track.RotateZ(-phi0_);
  // rotate about Y-axis to coordinate system that points along the central ray
  track.RotateY(-theta0_);
  return {asin(track.X() / part.momentum()), asin(track.Y() / part.momentum()),
          track.Z()};
}
ROOT::Math::PxPyPzMVector
spectrometer::detected_track(const particle& part, const double th_in,
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

void spectrometer::process(event& e) const {
  for (auto& part : e) {
    if (part.final_state() &&
        std::any_of(pid_.begin(), pid_.end(), [&](const auto& good_pid) {
          return part.type<int>() == good_pid;
        })) {
      LOG_JUNK2(name_, "Found matching final state particle " + part.name() +
                           " (status: " + std::to_string(part.status<int>()) +
                           ", momentum: " + std::to_string(part.momentum()) +
                           ")");
      // check spectrometer cuts
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
} // namespace liege
