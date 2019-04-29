#include "lp_gamma.hh"
#include <cmath>
#include <pcsim/core/particle.hh>
#include <pcsim/core/pdg.hh>

namespace pcsim {
namespace reconstruction {

namespace {
// TODO this helper function should really be part of configuration
bool get_bool(const configuration& conf, const string_path& path) {
  auto requirement = conf.get_optional<bool>(path);
  if (requirement && *requirement) {
    return true;
  }
  return false;
}
} // namespace

lp_gamma::lp_gamma(const configuration& conf, const string_path& path,
                   std::shared_ptr<TRandom> r)
    : lp_gamma::base_type{std::move(r)}
    , require_leading_{get_bool(conf, path / "require_leading")}
    , require_scat_{get_bool(conf, path / "require_scat")}
    , veto_leading_{get_bool(conf, path / "veto_leading")}
    , veto_scat_{get_bool(conf, path / "veto_scat")} {
  LOG_INFO("reconstruction",
           "Require leading hadron: " +
               std::string(require_leading_ ? "true" : "false"));
  LOG_INFO("reconstruction", "Require scattered lepton: " +
                                 std::string(require_scat_ ? "true" : "false"));
  LOG_INFO("reconstruction", "Veto leading hadron: " +
                                 std::string(veto_leading_ ? "true" : "false"));
  LOG_INFO("reconstruction", "Veto scattered lepton: " +
                                 std::string(veto_scat_ ? "true" : "false"));
  tassert(!(require_leading_ && veto_leading_),
          "Cannot require a reconstructed "
          "leading hadron but also veto the "
          "same particle");
  tassert(!(require_scat_ && veto_scat_), "Cannot require a reconstructed "
                                          "scattered lepton, but also veto the "
                                          "same particle");
}

void lp_gamma::process(lp_gamma_event& e) const {
  // do additional event reconstruction
  for (int i = 0; i < e.detected().size(); ++i) {
    // update the detected particle indices
    e.update_detected_index(i);
    // leading particle from decay products
    if (e.detected_leading_index() < 0 && e.leading_index() >= 0) {
      if (e.detected(i).generated().parent_first() == e.leading_index() &&
          e.detected(i).generated().n_parents() == 1 &&
          e.leading().n_daughters() > 1) {
        LOG_JUNK2("reconstruction",
                  "Found decay particle from leading particle, "
                  "scanning for the matching particles (particle type : " +
                      e.leading().name() + ") ");
        // disable enforced number of daughters for VMs, as RC particles are
        // also stored as daughters. For now we implicitly reconstruct using
        // only 2 particles for VMs
        if (e.leading().n_daughters() != 2 &&
            e.leading().type() != pdg_id::J_psi &&
            e.leading().type() != pdg_id::upsilon) {
          LOG_DEBUG("reconstruction",
                    "Only 2-body reconstruction implemented, but this "
                    "leading particle has " +
                        std::to_string(e.leading().n_daughters()) +
                        " daughters");
        } else {
          // locate the other decay product
          for (int j = i + 1; j < e.detected().size(); ++j) {
            if (e.detected(j).generated().parent_first() == e.leading_index()) {
              LOG_JUNK2("reconstruction",
                        "Found matching decay particle, performing "
                        "2-body reconstruction");
              // we found both!
              e.add_detected(
                  {e.leading(), e.detected(i).p() + e.detected(j).p(), -1});
              break;
            }
          }
        }
      }
    }
  }
  // check if we are fullfilling all requirments
  if (require_leading_ && e.detected_leading_index() < 0) {
    LOG_JUNK2("reconstruction",
              "Leading particle not reconstructed but required, "
              "setting event weight to zero.");
    e.update_weight(0);
  } else if (require_scat_ && e.detected_scat_index() < 0) {
    LOG_JUNK2("reconstruction",
              "Scattered lepton not reconstructed but required, "
              "setting event weight to zero.");
    e.update_weight(0);
  } else if (veto_leading_ && e.detected_leading_index() >= 0) {
    LOG_JUNK2("reconstruction", "Vetoing event with reconstructed leading "
                                "particle; setting event weight to zero.");
    e.update_weight(0);
  } else if (veto_scat_ && e.detected_scat_index() >= 0) {
    LOG_JUNK2("reconstruction", "Vetoing event with reconstructed scattered "
                                "lepton; setting event weight to zero.");
    e.update_weight(0);
  } else {
    LOG_JUNK2("reconstruction",
              "All reconstruction requirements fulfilled: good event!");
  }

  // that's all
}

} // namespace reconstruction
} // namespace pcsim
