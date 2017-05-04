#include "solid.hh"
#include <TMath.h>

namespace pcsim {
namespace detector {

#ifndef SOLID_ACCEPTANCE_PATH
#error SOLID_ACCEPTANCE_PATH NOT DEFINED
#endif

solid::solid(const configuration&, const string_path&,
             std::shared_ptr<TRandom> r)
    : solid::base_type{r} {
  acceptance_file_ = std::make_shared<TFile>(
      SOLID_ACCEPTANCE_PATH
      "acceptance_solid_JPsi_electron_target315_output.root",
      "read");
  tassert(acceptance_file_ && acceptance_file_->IsOpen(),
          "Error opening the acceptance file");
  forward_angle_ = static_cast<TH2F*>(
      acceptance_file_->Get("acceptance_ThetaP_forwardangle"));
  large_angle_ =
      static_cast<TH2F*>(acceptance_file_->Get("acceptance_ThetaP_largeangle"));
  tassert(forward_angle_ && large_angle_, "Failed to load acceptance histos");
}

void solid::process(event& e) const {
  for (auto& part : e) {
    if (part.final_state()) {
      LOG_JUNK2("solid",
                "Checking acceptance for final state particle " + part.name() +
                    " (status: " + std::to_string(part.status<int>()) + ")");
      const double theta = part.theta() * TMath::RadToDeg();
      const double p = part.momentum();
      LOG_JUNK2("solid", "theta [deg.]: " + std::to_string(theta) +
                             " p [GeV]: " + std::to_string(p));

      double acc = forward_angle_->GetBinContent(
          forward_angle_->GetXaxis()->FindBin(theta),
          forward_angle_->GetYaxis()->FindBin(p));
      int region = 0;

      if (acc > 0) {
        region = 1;
      } else {
        acc = large_angle_->GetBinContent(
            large_angle_->GetXaxis()->FindBin(theta),
            large_angle_->GetYaxis()->FindBin(p));
        if (acc > 0) {
          region = 2;
        }
      }
      LOG_JUNK2("solid", "Calculated acceptance: " + std::to_string(acc * 100) +
                             "\% (region: " + std::to_string(region) + ")");

      if (region == 0 || rng()->Uniform(0, 1) > acc) {
        LOG_JUNK2("solid", "Not accepted");
        // not accepted, continue
        continue;
      }
      LOG_JUNK2("solid", "Particle accepted!");

      e.add_detected({part, region});
    }
  }
}

} // namespace detector
} // namespace pcsim
