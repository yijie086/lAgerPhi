#ifndef PCSIM_BEAM_PHOTON_LOADED
#define PCSIM_BEAM_PHOTON_LOADED

#include <TRandom.h>
#include <cmath>
#include <memory>
#include <pcsim/core/factory.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/beam/beam.hh>
#include <pcsim/physics/flux.hh>

namespace pcsim {
namespace beam {

// =============================================================================
// beam::photon_data
//
// secondary photon beam
// =============================================================================
class photon_data  : data {
  // generalized constructors:
  // make collinear real photon event with energy E
  static photon_data make_real(const particle& lepton, const particle& target,
                               const double E, const double xs = 1.);

  // generate virtual photon event with Q2 and y
  // also needs access to a RNG to generate the azimuthal angle
  static photon_data make_virtual(const particle& lepton,
                                  const particle& target, const double Q2,
                                  const double y, std::shared_ptr<TRandom> rng,
                                  const double xs = 1.);

  photon_data() = default;
  photon_data(const photon_data&) = default;
  photon_data& operator=(const photon_data&) = default;

  photon_data(const double xs) : data{xs} {}

  photon_data(const particle::XYZTVector& p) : beam{{pdg_id::gamma, p}} {}
  photon_data(const particle::XYZTVector& p, const double xs)
      : beam{{pdg_id::gamma, p}, xs} {}

  double W2() const { return W2_; }
  double Q2() const { return Q2_; }
  double nu() const { return nu_; }
  double x() const { return x_; }
  double y() const { return y_; }

  const particle& scat() const { return scat_; }

private:
  double W2_{0.}; // invariant mass of photon-target system
  double Q2_{0.}; // photon virtuality
  double nu_{0.}; // photon enery in target rest frame
  double x_{0.};  // Bjorken x
  double y_{0.};  // energy fraction of photon in target rest frame
  particle scat_; // scattered lepton
};

// =============================================================================
// Base class for secondary photons from electron/positron beams on a nucleon
// Note: lepton beam is refered to as "beam", proton beam as "target"
// =============================================================================
class photon : public generator<photon_data, particle, particle> {
public:
  static factory<photon> factory;

  photon(std::shared_ptr<TRandom> r) : generator{std::move(r)} {}

  virtual photon_data generate(const beam::data& beam,
                               const beam::data& target) = 0;

private:
  ; // nothing here
};

// Bremsstrahlung photons
class bremsstrahlung : public photon {
public:
  enum class model { FLAT, PARAM, APPROX };

  bremsstrahlung(const configuration& cf, const string_path& path,
                 std::shared_ptr<TRandom> r);

  virtual photon_data generate(const beam::data& beam,
                               const beam::data& target);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return E_range_.width(); }

protected:
  double intensity(const double E, const double beam) const;

private:
  const model model_; // bremsstrahlung model
  const double rl_;   // number of radiation lenghts (when using approx model,
                      // set to zero otherwise)
  const double E_beam_;            // (maximum) electron beam energy
  const interval<double> E_range_; // photon energy range
  const double max_;               // the maximum intensity.
};

// virtual photons
class vphoton : public photon {
public:
  vphoton(const configuration& cf, const string_path& path,
          std::shared_ptr<TRandom> r);

  virtual photon_data generate(const beam::data& beam,
                               const beam::data& target);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const {
    return logy_range_.width() * logQ2_range_.width();
  }

protected:
  double flux(const double Q2, const double y, const particle& beam,
              const particle& target) const {
    return physics::flux::gamma_t_log(Q2, y, beam, target);
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


}
}

#endif
