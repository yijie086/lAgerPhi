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

#ifndef LAGER_PROC_DETECTOR_COMPOSITE_LOADED
#define LAGER_PROC_DETECTOR_COMPOSITE_LOADED

#include <memory>
#include <lager/core/configuration.hh>
#include <lager/proc/detector/detector.hh>
#include <vector>

namespace lager {
namespace detector {

class composite : public detector {
public:
  using base_type = detector;

  composite(const configuration&, const string_path&, std::shared_ptr<TRandom> r);

  virtual void process(event& e) const {
    for (const auto& det : detectors_) {
      det->process(e);
    }
  }

private:
  std::vector<std::shared_ptr<detector>> detectors_;
};

} // namespace detector
} // namespace lager

#endif
