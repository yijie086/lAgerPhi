#include "lp_gamma_generator.hh"
#include <pcsim/core/factory.hh>

namespace pcsim {

lp_gamma_generator::lp_gamma_generator(const configuration& cf,
                                       const string_path& path,
                                       std::shared_ptr<TRandom> r)
    : base_type{cf, path, r}
    , lepton_gen_{FACTORY_CREATE(beam::primary_generator, conf(), "beam", r)}
    , proton_gen_{FACTORY_CREATE(beam::primary_generator, conf(), "target", r)}
    , photon_gen_{FACTORY_CREATE(beam::photon_generator, conf(), "photon", r)}
    , decay_proc_{std::make_shared<decay::lp_gamma>(r)}
    , detector_proc_{FACTORY_CREATE(detector::detector, cf, "detector", r)}
    , rc_proc_{std::make_shared<reconstruction::lp_gamma>(r)} {
  register_initial(lepton_gen_);
  register_initial(proton_gen_);
  register_initial(photon_gen_);
}

lp_gamma_data lp_gamma_generator::generate_initial() const {
  auto beam = lepton_gen_->generate();
  auto target = proton_gen_->generate();
  auto photon = photon_gen_->generate(beam, target);
  if (photon.cross_section() <= 0) {
    return {0.};
  }
  return {beam, target, photon};
}

void lp_gamma_generator::build_event(lp_gamma_event& e) const {
  decay_proc_->process(e);
  detector_proc_->process(e);
  rc_proc_->process(e);
}

} // namespace pcsim
