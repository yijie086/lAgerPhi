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

#ifndef LAGER_CORE_SPLINE_LOADED
#define LAGER_CORE_SPLINE_LOADED

#include <TSpline.h>
#include <lager/core/configuration.hh>
#include <lager/core/interval.hh>

namespace lager {

// a 1D spline interpolation that uses datapoints obtained from a configuration
// object
class spline : public configurable {
public:
  spline(const configuration& conf, const string_path& path)
      : configurable{conf, path}
      , x_{conf().get_vector<double>("x")}
      , y_{conf().get_vector<double>("y")}
      , spline_{path.str().c_str(), &x_[0], &y_[0], static_cast<int>(x_.size())}
      , xrange_{x_.front(), x_.back()} {}

  double eval(const double x) const {
    if (xrange_.includes(x)) {
      return spline_.Eval(x);
    } else {
      return 0.;
    }
  }
  double operator()(const double x) const { return eval(x); }

private:
  std::vector<double> x_;
  std::vector<double> y_;
  TSpline3 spline_;
  interval<double> xrange_;
};
}

#endif
