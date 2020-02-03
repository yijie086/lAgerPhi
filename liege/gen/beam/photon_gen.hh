#ifndef LIEGE_GEN_BEAM_PHOTON_GEN_LOADED
#define LIEGE_GEN_BEAM_PHOTON_GEN_LOADED

#include <TRandom.h>
#include <memory>
#include <liege/core/factory.hh>
#include <liege/core/generator.hh>
#include <liege/core/particle.hh>
#include <liege/gen/beam/generator.hh>
#include <liege/gen/beam/photon.hh>
#include <liege/gen/beam/primary.hh>
#include <liege/physics/photon.hh>

namespace liege {
namespace beam {

// Bremsstrahlung photons
class bremsstrahlung : public photon_generator {
public:
  enum class model { FLAT, PARAM, APPROX, EXACT };

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

// Bremsstrahlung photons for a realistic (extended) target
class bremsstrahlung_realistic_target : public photon_generator {
public:
  bremsstrahlung_realistic_target(const configuration& cf,
                                  const string_path& path,
                                  std::shared_ptr<TRandom> r);

  virtual photon generate(const primary& beam, const primary& target);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return E_range_.width(); }

protected:
  double intensity(const double E, const double beam, const double vz) const;
  double total_rl(const double vz) const;

private:
  const double rl_radiator_;            // RL for the radiator
  const double rl_window_;              // RL for the target window
  const double extra_rl_per_cm_;        // extra RL/cm in the target
  const interval<double> target_range_; // target range in cm

  const double E_beam_;            // (maximum) electron beam energy
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
} // namespace liege

#endif
