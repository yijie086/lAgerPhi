#ifndef PCSIM_GEN_PHOTON_LOADED
#define PCSIM_GEN_PHOTON_LOADED

#include <TRandom.h>
#include <cmath>
#include <memory>
#include <pcsim/core/factory.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/beam.hh>
#include <pcsim/physics/flux.hh>

namespace pcsim {
namespace gen {

struct photon_data  : generator_data {
  double W2;       // invariant mass of photon-target system
  double Q2;       // photon virtuality
  double nu;       //
  double x;
  double y;
  particle scat;   // scattered lepton
  particle photon; // photon (virtual or real)
  photon_data() = default;
  photon_data(double xs, double ps) : generator_data(xs, ps) {}
};

// =============================================================================
// Base class for secondary photons from electron/positron beams on a nucleon
// Note: lepton beam is refered to as "beam", proton beam as "target"
// =============================================================================
class photon : public generator<photon_data, particle, particle> {
public:
  static factory<photon> factory;

  photon(const configuration& conf, const string_path& path,
         std::shared_ptr<TRandom> r)
      : generator{conf, path, std::move(r)} {}

  virtual photon_data generate(const particle& beam,
                               const particle& target) = 0;

private:
  ; // nothing here
};

// Bremsstrahlung photons
class bremsstrahlung : public photon {
public:
  enum class model { FLAT, PARAM, APPROX };

  bremsstrahlung(const configuration& conf, const string_path& path,
                 std::shared_ptr<TRandom> r);

  virtual photon_data generate(const particle& beam, const particle& target);
  virtual double max_cross_section() const { return max_; }

protected:
  double intensity(const double E, const double beam) const;

private:
  const model model_; // bremsstrahlung model
  const double rl_;   // number of radiation lenghts (when using approx model,
                      // set to zero otherwise)
  const double E_beam_; // (maximum) electron beam energy
  const interval<double> E_range_; // photon energy range
  const double max_;               // the maximum intensity.
};

// virtual photons
class vphoton : public photon {
public:
  vphoton(const configuration& cf, const string_path& path,
          std::shared_ptr<TRandom> r);

  virtual photon_data generate(const particle& beam, const particle& target);
  virtual double max_cross_section() const { return max_; }

protected:
  double flux(const double Q2, const double y, const particle& beam,
              const particle& target) const {
    return physics::flux::gamma_t_log(Q2, y, beam.mom, target.mom);
  }

private:
  interval<double> Q2_range(const particle& beam, const particle& target,
                            const double y) const;
  double calc_max_flux(const configuration& cf) const;
  interval<double> calc_max_Q2_range(const configuration& cf) const;
  particle generate_scat(const double Q2, const double y, particle beam,
                         const particle& target);

  // primary kinematic boundaries
  const interval<double> y_range_;
  const interval<double> Q2_range_;
  // derived kinematic boundaries
  const interval<double> logy_range_;
  const interval<double> logQ2_range_;

  // maximum flux
  const double max_;
};


}
}

#endif
