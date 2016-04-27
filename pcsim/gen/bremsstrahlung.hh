#ifndef PCSIM_GEN_BREMSSTRAHLUNG_LOADED
#define PCSIM_GEN_BREMSSTRAHLUNG_LOADED

#include <TF1.h>
#include <pcsim/core/accept_reject.hh>
#include <pcsim/core/configuration.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/interval.hh>
#include <pcsim/physics/bremsstrahlung.hh>

namespace pcsim {
namespace gen {

struct photon_beam {
  double energy; // Photon energy
  double flux;   // additional factor to relate bremsstrahlung distribution to
                 // incoming electron flux
  bool good;     // Is this a good photon?
  constexpr photon_beam(double E, double f, bool g = true)
      : energy{E}, flux{f}, good{g} {}
};

class bremsstrahlung : public generator<bremsstrahlung, photon_beam> {
public:
  using base_type = generator<bremsstrahlung, photon_beam>;

  bremsstrahlung(const ptree& settings, const string_path& path,
                 std::shared_ptr<TRandom> r)
      : base_type{settings, path, "Bremsstrahlung Generator",
                  std::move(r)}
      , E0_{conf().get<double>("electron_energy")}
      , range_{conf().get_range<double>("range")}
      , integral_{calc_integral()}
      , imax_{intensity(range_.min)} {}

  photon_beam gen_impl() {
    double E = accept_reject_1D(
        rng(), range_, [=](const double k) { return intensity(k); }, imax_);
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
