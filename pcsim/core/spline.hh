#ifndef PCSIM_CORE_SPLINE_LOADED
#define PCSIM_CORE_SPLINE_LOADED

#include <TSpline.h>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/interval.hh>

namespace pcsim {

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
