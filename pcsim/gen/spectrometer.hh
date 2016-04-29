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

namespace pcsim {
namespace gen {

struct spec_track {
  bool accept;
};

class spectrometer : public generator<spectrometer, spec_track> {
public:
  using base_type = generator<spectrometer, spec_track>;

  spectrometer(const ptree& settings, const string_path& path,
               std::shared_ptr<TRandom> r);

  spec_track check(TLorentzVector track) {
    // rotate the track to the spectrometer center
    track.RotateY(-theta_);

    const double p = track.Vect().Mag();
    // get the angles
    const double sx = track.X() / p;
    const double sy = track.Y() / p;
    const double thx = std::fabs(asin(sx));
    const double thy = std::fabs(asin(sy));

    return {p_range_.includes(p) && thx < x_acc_ && thy < y_acc_};
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
