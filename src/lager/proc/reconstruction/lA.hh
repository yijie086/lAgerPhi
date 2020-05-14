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

#ifndef LAGER_PROC_RECONSTRUCTION_LA_LOADED
#define LAGER_PROC_RECONSTRUCTION_LA_LOADED

#include <lager/gen/lA_event.hh>
#include <lager/proc/reconstruction/reconstruction.hh>

namespace lager {
namespace reconstruction {
// =============================================================================
// RECONSTRUCT ADDITIONAL PARTICLES IN A GAMMA_P EVENT
//
// also handles trigger-level cuts
//
// Note: sets the event weight to zero for events that don't pass the cut
//       (will then be removed by the event_generator)
// =============================================================================

class lA : public reconstruction<lA_event> {
public:
  using base_type = reconstruction<lA_event>;
  lA(const configuration& conf, const string_path& path,
           std::shared_ptr<TRandom> r);
  virtual void process(lA_event& e) const;

private:
  const bool require_leading_{false};
  const bool veto_leading_{false};
  const bool require_scat_{false};
  const bool veto_scat_{false};
  const bool require_recoil_{false};
  const bool veto_recoil_{false};
};

} // namespace reconstruction
} // namespace lager

#endif
