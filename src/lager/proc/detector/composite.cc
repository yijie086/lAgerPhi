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

#include "composite.hh"
#include <lager/core/factory.hh>

namespace lager {
namespace detector {

composite::composite(const configuration& cf, const string_path& path,
                     std::shared_ptr<TRandom> r)
    : composite::base_type{r} {
  auto tmp_conf = cf;
  auto& conf = tmp_conf.raw_node(path / "components");
  for (const auto& child : conf) {
    string_path child_path = path / "components" / child.first.c_str();
    LOG_INFO("composite", "Constructing child detector: " + child_path.str());
    auto det = FACTORY_CREATE(detector, tmp_conf, child_path, r);
    detectors_.push_back(det);
  }
}
} // namespace detector

} // namespace lager
