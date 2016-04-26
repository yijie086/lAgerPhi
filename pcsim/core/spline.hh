#ifndef PCSIM_CORE_SPLINE_LOADED
#define PCSIM_CORE_SPLINE_LOADED

#include <Spline3.h>
#include <pcsim/util/configuration.hh>
#include <pcsim/util/interval.hh>

namespace pcsim {

// a 1D spline interpolation that uses datapoints obtained from a configuration
// object
class spline : public configurable {
public:
  spline(const ptree& settings, const string_path& path)
      : configurable{settings, path}
      , x_{conf().get_vector<double>(path + "/x")}
      , y_{conf().get_vector<double>(path + "/y")}
      , spline_{path.str().c_str(), &x_[0], &y_[0], x_.size()}
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
  const std::vector<double> x_;
  const std::vector<double> y_;
  const TSpline3 spline_;
  const interval<double> xrange_;
};
}

#endif
