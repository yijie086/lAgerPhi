#include "fixed_target.hh"

namespace lager::initial {

realistic_target::realistic_target(const configuration& cf,
                                   const string_path& path)
    : rl_pre_{calc_rl_pre(cf, path)}
    , rl_per_cm_{1.0 / (cf.get<double>(path / "target" / "radiation_length") /
                        cf.get<double>(path / "target" / "density"))}
    , target_range_{cf.get_range<double>(path / "target" / "range")} {
  // radiation length
  LOG_INFO("realistic_target",
           "Total pre (rad+window) RL [%]: " + std::to_string(rl_pre_ * 100.));
  LOG_INFO("realistic_target", "Target RL (for 1cm of target material) [%]: " +
                                   std::to_string(rl_per_cm_ * 100.));
  LOG_INFO("realistic_target",
           "Target RL (for all target material) [%]: " +
               std::to_string(rl_per_cm_ * target_range_.width() * 100.));
  LOG_INFO("realistic_target",
           "RL at front of target [%]" +
               std::to_string(total_rl(target_range_.min) * 100.));
  LOG_INFO("realistic_target",
           "RL at back of target [%]" +
               std::to_string(total_rl(target_range_.max) * 100.));
}

double realistic_target::calc_rl_pre(const configuration& cf,
                                     const string_path& path) const {
  const double rl_rad =
      cf.get<double>(path / "radiator" / "thickness") /
      (cf.get<double>(path / "radiator" / "radiation_length") /
       cf.get<double>(path / "radiator" / "density"));
  const double rl_window =
      cf.get<double>(path / "window" / "thickness") /
      (cf.get<double>(path / "window" / "radiation_length") /
       cf.get<double>(path / "window" / "density"));
  LOG_INFO("realistic_target",
           "Radiator RL (calculated) [%]: " + std::to_string(rl_rad * 100.));
  LOG_INFO("realistic_target",
           "Window RL (calculated) [%]: " + std::to_string(rl_window * 100.));
  return rl_rad + rl_window;
}
} // namespace lager::initial
