// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
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

#ifndef PCSIM_GEN_LA_GENERATOR_LOADED
#define PCSIM_GEN_LA_GENERATOR_LOADED

#include <lager/core/generator.hh>
#include <lager/gen/initial/generator.hh>
#include <lager/gen/lA/generator.hh>
#include <lager/gen/lA_event.hh>
#include <lager/proc/decay/lA.hh>
#include <lager/proc/detector/detector.hh>
#include <lager/proc/reconstruction/lA.hh>
#include <memory>

namespace lager {

class lA_generator : public event_generator<lA_event, lA_data> {
public:
  using base_type = event_generator<lA_event, lA_data>;

  lA_generator(const configuration& cf, const string_path& path,
               std::shared_ptr<TRandom> r);

protected:
  virtual lA_data generate_initial() const;
  virtual void build_event(lA_event& e) const;

private:
  // initial state generators
  std::shared_ptr<initial::vertex_generator> vertex_gen_;
  std::shared_ptr<initial::beam_generator> lepton_gen_;
  std::shared_ptr<initial::beam_generator> ion_gen_;
  std::shared_ptr<initial::target_generator> target_gen_;
  std::shared_ptr<initial::photon_generator> photon_gen_;
  // event processors
  std::shared_ptr<decay::lA> decay_proc_;
  std::shared_ptr<detector::detector> detector_proc_;
  std::shared_ptr<reconstruction::lA> rc_proc_;
};

} // namespace lager

#endif
