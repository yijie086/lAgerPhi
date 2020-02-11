// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2020 Sylvester Joosten <sjoosten@anl.gov>
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

#include "lA_generator.hh"
#include <lager/core/factory.hh>

namespace lager {

lA_generator::lA_generator(const configuration& cf, const string_path& path,
                           std::shared_ptr<TRandom> r)
    : base_type{cf, path, r}
    , vertex_gen_{FACTORY_CREATE(initial::vertex_generator, conf(), "vertex",
                                 r)}
    , lepton_gen_{FACTORY_CREATE(initial::beam_generator, conf(), "beam/lepton",
                                 r)}
    , ion_gen_{FACTORY_CREATE(initial::beam_generator, conf(), "beam/ion", r)}
    , target_gen_{FACTORY_CREATE(initial::target_generator, conf(), "target",
                                 r)}
    , photon_gen_{FACTORY_CREATE(initial::photon_generator, conf(), "photon",
                                 r)}
    , decay_proc_{std::make_shared<decay::lA>(cf, "decay", r)}
    , detector_proc_{FACTORY_CREATE(detector::detector, cf, "detector", r)}
    , rc_proc_{std::make_shared<reconstruction::lA>(cf, "reconstruction", r)} {
  register_initial(lepton_gen_);
  register_initial(ion_gen_);
  register_initial(target_gen_);
  register_initial(photon_gen_);
}

lA_data lA_generator::generate_initial() const {
  auto vertex = vertex_gen_->generate();
  auto lepton = lepton_gen_->generate(vertex);
  auto ion = ion_gen_->generate(vertex);
  auto target = target_gen_->generate(ion);
  if (target.cross_section() <= 0) {
    return {0.};
  }
  auto photon = photon_gen_->generate(lepton, target);
  if (photon.cross_section() <= 0) {
    return {0.};
  }
  return {lepton, ion, target, photon};
}

void lA_generator::build_event(lA_event& e) const {
  decay_proc_->process(e);
  detector_proc_->process(e);
  rc_proc_->process(e);
}

} // namespace lager
