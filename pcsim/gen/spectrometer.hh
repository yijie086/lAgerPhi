#ifndef PCSIM_GEN_SPECTROMETER_LOADED
#define PCSIM_GEN_SPECTROMETER_LOADED

#include <TLorentzVector.h>
#include <TRandom.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/interval.hh>
#include <pcsim/physics/jpsi.hh>
#include <string>

#include <pcsim/core/assert.hh>
#include <pcsim/core/stringify.hh>

namespace pcsim {
namespace gen {

struct spec_track {
  TLorentzVector track;
  int charge = 0;
  double p = 0;
  double thx = 0;
  double thy = 0;
  bool accept = false;

  spec_track() {}
  spec_track(const TLorentzVector& track, const int charge)
      : track{track}
      , charge{charge}
      , p{track.Vect().Mag()}
      , thx{asin(track.X() / p)}
      , thy{asin(track.Y() / p)} {}
};

class spectrometer : public generator<spectrometer, spec_track> {
public:
  using base_type = generator<spectrometer, spec_track>;

  spectrometer(const configuration& conf, const string_path& path,
               std::shared_ptr<TRandom> r);

  spec_track check(TLorentzVector track, const int charge) {

    // rotate to spectrometer system
    track.RotateY(-theta_);
    // create spectrometer track
    spec_track t{track, charge};
    // do we accept this track?
    t.accept = p_range_.includes(t.p) && std::fabs(t.thx) < x_acc_ &&
               std::fabs(t.thy) < y_acc_;
    // that's all!

    static int debug = 0;
    if (false && t.p > 6 && p_range_.includes(t.p)) {
      std::cout << "\n"
                << t.p << " " << p_range_.min << " " << p_range_.max << " ("
                << t.thx << "," << t.thy << ") (" << x_acc_ << "," << y_acc_
                << ") " << (t.accept ? "accept" : "reject") << "\n";
      tassert(++debug < 100, "DEBUG");
    }
    return t;
  }


private:
  const interval<double> p_range_;
  const double theta_;
  const double x_acc_;
  const double y_acc_;

};

} // gen
} // pcsim

#endif
