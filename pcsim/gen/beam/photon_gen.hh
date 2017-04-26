#ifndef PCSIM_GEN_BEAM_PHOTON_GEN_LOADED
#define PCSIM_GEN_BEAM_PHOTON_GEN_LOADED

#include <TRandom.h>
#include <memory>
#include <pcsim/core/factory.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/beam/beam.hh>
#include <pcsim/gen/beam/photon.hh>
#include <pcsim/gen/beam/primary.hh>
#include <pcsim/physics/photon.hh>

namespace pcsim {
namespace beam {

// Bremsstrahlung photons
class bremsstrahlung : public photon_generator {
public:
  enum class model { FLAT, PARAM, APPROX };

  bremsstrahlung(const configuration& cf, const string_path& path,
                 std::shared_ptr<TRandom> r);

  virtual photon generate(const primary& beam, const primary& target);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return E_range_.width(); }

protected:
  double intensity(const double E, const double beam) const;

private:
  const model model_;   // bremsstrahlung model
  const double rl_;     // number of radiation lenghts (when using approx model,
                        // set to zero otherwise)
  const double E_beam_; // (maximum) electron beam energy
  const interval<double> E_range_; // photon energy range
  const double max_;               // the maximum intensity.
};

// virtual photons
class vphoton : public photon_generator {
public:
  vphoton(const configuration& cf, const string_path& path,
          std::shared_ptr<TRandom> r);

  virtual photon generate(const primary& beam, const primary& target);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const {
    return logy_range_.width() * logQ2_range_.width();
  }

protected:
  double flux(const double Q2, const double y, const particle& beam,
              const particle& target) const {
    return physics::gamma_t_log(Q2, y, beam, target);
  }

private:
  double calc_max_flux(const configuration& cf) const;
  interval<double> calc_max_Q2_range(const configuration& cf) const;
  interval<double> calc_max_W2_range(const configuration& cf) const;

  // primary kinematic boundaries
  const interval<double> y_range_;
  const interval<double> Q2_range_;
  // derived kinematic boundaries
  const interval<double> logy_range_;
  const interval<double> logQ2_range_;
  // additional cuts
  const interval<double> W2_range_;

  // maximum flux
  const double max_;
};

} // namespace beam
} // namespace pcsim

#endif
