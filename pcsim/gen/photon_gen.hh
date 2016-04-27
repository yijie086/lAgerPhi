#ifndef PCSIM_PHOTON_GEN_LOADED
#define PCSIM_PHOTON_GEN_LOADED

#include <TF1.h>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/interval.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/physics/photon_beam.hh>

namespace pcsim {

struct photon_beam {
  double energy; // Photon energy
  double weight; // additional event weight
  bool good;     // Is this a good photon?
  constexpr photon_beam(double E, double w, bool g = true)
      : energy{E}, weight{w}, good{g} {}
}


class photon_gen : public generator<photon_gen> {
public:
  using event_type = photon_beam;
  using base_type = generator<photon_gen>;

  photon_gen(const ptree& settings, const string_path& path,
             std::shared_ptr<TRandom> r)
      : base_type{settings, path, "Secondary Photon Beam Generator",
                  std::move{r}}
      , E0_{conf().get_range<double>("electron_energy")}
      , range_{conf().get_range<double>("range")}
      , integral_{calc_integral()}
      , imax_{intensity(range_.min)} {}

  photon_beam gen_impl() {
    // get an energy candidate
    double E = rng().Uniform(range_.min, range_.max);
    // check if we should accept or reject this number
    double test = rng().Uniform(0, imax_);
    if (test > intensity(E)) {
      // reject --> try again
      return gen_impl();
    }
    return {E, integral_};
  }

private:
  double intensity(const double k) {
    return photon_intensity_10_param(E0_, x[0]);
  }
  double calc_integral() {
    auto f = [](const double* x, const double* p) { return intensity(x[0]); };
    TF1 tf("brems10", f, range_.min, range_.max);
    return tf.Integral(range_.min, range_.max);
  }

  const double E0_;              // primary electron beam energy
  const interval<double> range_; // photon energy range
  const double integral_;        // integrated photon intensity within the range
  const double imax_; // the maximum intensity. Bremsstrahlung spectrum is
                      // monotonously falling, and therefor imax is given by
                      // f(range_.min)
}



}

#endif
