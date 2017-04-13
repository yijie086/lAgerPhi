#ifndef PCSIM_GEN_PROCESS_GAMMA_P_1X_LOADED
#define PCSIM_GEN_PROCESS_GAMMA_P_1X_LOADED
#define PCSIM_GEN_PROCESS_GAMMA_P_LOADED

#include <pcsim/core/factory.hh>
#include <pcsim/core/pdg.hh>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/beam/beam.hh>
#include <pcsim/gen/beam/photon.hh>

namespace pcsim {
namespace process {

// =============================================================================
// gamma_p_1X_data
//
// resonant production data
// =============================================================================
class gamma_p_1X_data : public generator_data {
public:
  gamma_p_1X_data() = default;
  gamma_p_1X_data(const gamma_p_1X_data&) = default;
  gamma_p_1X_data& operator=(const gamma_p_1X_data&) = default;

  gamma_p_1X_data(const double xs) : generator_data{xs} {}

  gamma_p_1X_data(const beam::photon_data& photon, const beam::data& target,
                  const particle& X1, const double xs, const double R = 0);

  double R() const { return R_; }
  const particle& X() { return X_; };

private:
  double R_;
  particle X_;
};

// =============================================================================
// process::gamma_p_1X
//
// abstract base class for all gamma_p_1X processes
// =============================================================================
class gamma_p_1X
    : public generator<gamma_p_1X_data, beam::photon_data, beam::data> {
public:
  static factory<gamma_p_1X> factory;

  gamma_p_1X(std::shared_ptr<TRandom> r) : generator{std::move(r)} {}
};

// =============================================================================
// process::gamma_p_1X_qpq
//
// gamma + p -> Pq process
//
// (resontant production of a quarkonium pentaquark)
//
// Uses the following expressions (cf. pcsim/physics/vm.hh)
// (assuming the photo-production happens through a given quarkonium pole)
//  * R (sigma_L/sigma_T):
//        R_martynov(...)
//  * Dipole FF for sigma_gamma -> sigma_t:
//        dipole_ff__hermes(...)
//  * photo-production cross section:
//        simple gaussian
// =============================================================================
class gamma_p_1X_qpq : public gamma_p_1X {
public:
  gamma_p_1X_qpq(const configuration& cf, const string_path& path,
                       std::shared_ptr<TRandom> r);
  virtual gamma_p_1X_data generate(const beam::photon_data& photon,
                                   const beam::data& target);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return 1.; }

private:
  double calc_max_xsec(const configuration& cf) const;

  // cross section component evaluation
  double sigma(const double W2) const;
  double R(const double Q2) const;
  double dipole(const double Q2) const;

  // recoil and  particle info
  const particle vm_pole_;
  const particle qpq_;

  // cross section settings
  const double qpq_amplitude_; // peak cross section amplitude
  const double qpq_coupling_;  // coupling through the quarkonium pole
  const double R_vm_c_;        // c-parameter for R
  const double R_vm_n_;        // n-parameter for R
  const double dipole_n_;      // n-parameter for dipole factor

  const double max_;
};

} // namespace process
} // namespace pcsim

#endif
