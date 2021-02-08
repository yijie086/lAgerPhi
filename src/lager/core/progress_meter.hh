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

#ifndef LAGER_CORE_PROGRESS_METER_LOADED
#define LAGER_CORE_PROGRESS_METER_LOADED

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>

// =============================================================================
// A simple console progress meter
// =============================================================================
namespace lager {
class progress_meter {
public:
  constexpr static const size_t PRECISION = 10000.; // 0.01 percent
  explicit progress_meter(size_t max, size_t start_index = 0,
                          size_t precision = PRECISION)
      : max_{max}
      , index_{start_index}
      , precision_{precision > max ? max : precision} {
    std::cerr << "\nProcessing... " << std::endl;
    update();
  }
  void update(size_t i) {
    index_ = i;
    if (index_ > max_) {
      index_ = max_;
    }
    // update when needed
    // commented out because of FPE, need to fix this TODO
//    if (!index_ || !(index_ % (max_ / precision_) || !(index_ % 1000))) {
    if (!(index_ % 100)) {
      const double cnt = (index_ > 50) ? (100. * index_) / max_ : 0.;
      char msg[15];
      sprintf(msg, "  %3.2f%%\r", cnt);
      std::cerr << msg << std::flush;
    }
  }
  void update(size_t i, const size_t max) {
    max_ = max;
    update(i);
  }
  void update() {
    index_ += 1;
    update(index_);
  }
  ~progress_meter() { std::cerr << "      Done!" << std::endl; }

private:
  size_t max_;
  size_t index_;
  const size_t precision_;
};
} // ns lager

#endif
