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

#ifndef LAGER_PROC_DECAY_LA_LOADED
#define LAGER_PROC_DECAY_LA_LOADED

#include <lager/core/particle.hh>
#include <lager/gen/lA_event.hh>
#include <lager/physics/decay.hh>
#include <lager/proc/decay/decay.hh>

namespace lager {
namespace decay {
// =============================================================================
// DECAY ALL UNSTABLE PARTICLES IN A GAMMA_P EVENT
//
// also handles chained decays
//
// Supported channels:
//  * Pc according to Wang
//  * e+e- decay of VMs
//
// Note: adds event weight to account for branching ratios when not simulating
// the full decay width
// =============================================================================
class radiative_decay_vm {
public:
  radiative_decay_vm();
  void process(lA_event& e, const int vm_index);
};

class lA : public decay<lA_event> {
public:
  using base_type = decay<lA_event>;
  lA(const configuration&, const string_path&, std::shared_ptr<TRandom> r);
  virtual void process(lA_event& e) const;

private:
  void quarkonium_schc(lA_event& e, const int index) const;
  void pentaquark_qpq(lA_event& e, const int index) const;

  const particle vm_decay_lplus_;
  const particle vm_decay_lminus_;
  const double vm_decay_br_;
  std::unique_ptr<radiative_decay_vm> radiative_decay_;
};

} // namespace decay
} // namespace lager

#endif
