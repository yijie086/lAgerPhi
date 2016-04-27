#ifndef PCSIM_GENERATOR_BREMSSTRAHLUNG_LOADED
#define PCSIM_GENERATOR_BREMSSTRAHLUNG_LOADED

#include <TF1.h>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/interval.hh>
#include <pcsim/physics/bremsstrahlung.hh>

namespace pcsim {
namespace gen {

struct photon_beam {
  double energy; // Photon energy
  double weight; // additional event weight
  bool good;     // Is this a good photon?
  constexpr photon_beam(double E, double w, bool g = true)
      : energy{E}, weight{w}, good{g} {}
};

class bremsstrahlung : public generator<bremsstrahlung, photon_beam> {
public:
  using base_type = generator<bremsstrahlung, photon_beam>;

  bremsstrahlung(const ptree& settings, const string_path& path,
                 std::shared_ptr<TRandom> r)
      : base_type{settings, path, "Secondary Photon Beam Generator",
                  std::move(r)}
      , E0_{conf().get<double>("electron_energy")}
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

  const interval<double>& range() const { return range_; }

private:
  double intensity(const double k) {
    return physics::bremsstrahlung_intensity_10_param(E0_, k);
  }
  double calc_integral() {
    auto f = [=](double* x, double* p) {
      return physics::bremsstrahlung_intensity_10_param(E0_, x[0]);
    };
    TF1 tf("brems10", f, range_.min, range_.max, 0);
    return tf.Integral(range_.min, range_.max);
  }

  const double E0_;              // primary electron beam energy
  const interval<double> range_; // photon energy range
  const double integral_;        // integrated photon intensity within the range
  const double imax_; // the maximum intensity. Bremsstrahlung spectrum is
                      // monotonously falling, and therefor imax is given by
                      // f(range_.min)
};

} // gen
} // pcsim

#endif
