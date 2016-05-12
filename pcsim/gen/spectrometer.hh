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

    // apply smearing if we accepted the track (avoiding unnecessary work)
    if (t.accept) {
      // smear kinematics
      if (p_smear_>0) {
        t.p = rng()->Gaus(t.p, p_smear_ * t.p);
      }
      if (x_smear_ > 0) {
        t.thx = rng()->Gaus(t.thx, x_smear_);
      }
      if (y_smear_ > 0) {
        t.thy = rng()->Gaus(t.thy, y_smear_);
      }
      // recalculate track
      const double px = t.p * sin(t.thx);
      const double py = t.p * sin(t.thy);
      const double pz = sqrt(t.p * t.p - px * px - py * py);
      t.track.SetXYZM(px, py, pz, t.track.M());
      // rotate back to lab frame
      track.RotateY(theta_);
    }

    // that's all!

    return t;
  }


private:
  const interval<double> p_range_;
  const double theta_;
  const int charge_;
  const double x_acc_;
  const double y_acc_;
  const double p_smear_;
  const double x_smear_;
  const double y_smear_;
};

} // gen
} // pcsim

#endif
