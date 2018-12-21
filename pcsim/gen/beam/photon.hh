#ifndef PCSIM_GEN_BEAM_PHOTON_LOADED
#define PCSIM_GEN_BEAM_PHOTON_LOADED

#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/beam/primary.hh>

namespace pcsim {
namespace beam {

// =============================================================================
// beam::photon
//
// secondary photon beam
// =============================================================================
class photon : public primary {
public:
  // generalized constructors:
  // make collinear real photon event with energy E
  static photon make_real(const particle& lepton, const particle& target,
                          const double E, const double xs);

  // generate virtual photon event with Q2 and y
  // also needs an azimuthal angle in the CM frame
  static photon make_virtual(const particle& lepton, const particle& target,
                             const double Q2, const double y, const double xs,
                             const double phi);

  photon() = default;
  photon(const photon&) = default;
  photon& operator=(const photon&) = default;

  photon(const double xs) : primary{xs} {}

  photon(const particle::XYZTVector& p)
      : primary{{pdg_id::gamma, p, particle::status_code::SECONDARY_BEAM}} {}
  photon(const particle::XYZTVector& p, const double xs)
      : primary{{pdg_id::gamma, p, particle::status_code::SECONDARY_BEAM}, xs} {
  }

  double epsilon() const { return epsilon_; }
  double W() const { return sqrt(W2_); }
  double W2() const { return W2_; }
  double Q2() const { return Q2_; }
  double nu() const { return nu_; }
  double x() const { return x_; }
  double y() const { return y_; }

  const particle& scat() const { return scat_; }

private:
  double epsilon_{0.}; // gamma_L/gamma_T (xs stores just gamma_T)
  double W2_{0.};      // invariant mass of photon-target system
  double Q2_{0.};      // photon virtuality
  double nu_{0.};      // photon enery in target rest frame
  double x_{0.};       // Bjorken x
  double y_{0.};       // energy fraction of photon in target rest frame
  particle scat_;      // scattered lepton
};

} // namespace beam
} // namespace pcsim

#endif
