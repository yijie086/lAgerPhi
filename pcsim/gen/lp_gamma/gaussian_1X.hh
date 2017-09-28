#ifndef PCSIM_GEN_LP_GAMMA_GAUSSIAN_1X_LOADED
#define PCSIM_GEN_LP_GAMMA_GAUSSIAN_1X_LOADED

#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/lp_gamma/generator.hh>
#include <pcsim/gen/lp_gamma_event.hh>

namespace pcsim {
namespace lp_gamma {

// =============================================================================
// lp_gamma::gaussian_qpq
//
// gamma + p -> Pq process
//
// (resonant production of a quarkonium pentaquark)
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
class gaussian_qpq : public lp_gamma::generator {
public:
  using base_type = lp_gamma::generator;
  gaussian_qpq(const configuration& cf, const string_path& path,
               std::shared_ptr<TRandom> r);
  virtual lp_gamma_event generate(const lp_gamma_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return 1.; }

private:
  double calc_max_xsec(const configuration& cf) const;
  interval<double> calc_W2_range(const double n_sigma) const;

  // cross section component evaluation
  double sigma(const double W2) const;
  double R(const double Q2) const;
  double dipole(const double Q2) const;

  // recoil and  particle info
  const particle vm_pole_;          // relevant VM pole
  const particle qpq_;              // quarkonium pentaquark assumption
  const interval<double> W2_range_; // mass squared range around the pole mass

  // cross section settings
  const double qpq_amplitude_; // peak cross section amplitude
  const double qpq_coupling_;  // coupling through the quarkonium pole
  const double R_vm_c_;        // c-parameter for R
  const double R_vm_n_;        // n-parameter for R
  const double dipole_n_;      // n-parameter for dipole factor

  const double max_;
};

} // namespace lp_gamma
} // namespace pcsim

#endif
