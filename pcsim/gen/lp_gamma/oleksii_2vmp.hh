#ifndef PCSIM_GEN_LP_GAMMA_OLEKSII_2VMp_LOADED
#define PCSIM_GEN_LP_GAMMA_OLEKSII_2VMp_LOADED

#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/lp_gamma/generator.hh>
#include <pcsim/gen/lp_gamma_event.hh>

namespace pcsim {
namespace lp_gamma {

// =============================================================================
// lp_gamma::oleksii_2vmp
//
// gamma + p -> VM + p process
//
// Uses the following expressions (cf. pcsim/physics/vm.hh)
//  * R (sigma_L/sigma_T):
//        R_vm_martynov(...)
//  * Dipole FF for sigma_gamma -> sigma_t:
//        dipole_ff_vm_hermes(...)
//  * t-channel cross section:
//        Oleksii's implementation for J/psi and Upsilon
// =============================================================================

class oleksii_2vmp_amplitude;
class oleksii_2vmp_slope;

class oleksii_2vmp : public lp_gamma::generator {
public:
  using base_type = lp_gamma::generator;

  oleksii_2vmp(const configuration& cf, const string_path& path,
               std::shared_ptr<TRandom> r);
  virtual lp_gamma_event generate(const lp_gamma_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return max_t_range_.width(); }

private:
  double calc_max_b(const configuration& cf) const;
  double calc_max_xsec(const configuration& cf) const;
  interval<double> calc_max_t_range(const configuration& cf) const;

  // kinematics range
  interval<double> t_range(const double W2, const double Q2,
                           const double Mt) const;

  // cross section component evaluation
  double dsigma_dt(const double W2, const double t, const double b) const;
  double R(const double Q2) const;
  double dipole(const double Q2) const;

  // jacobian (equal to unity)
  double jacobian() const;

  // threshold squared for these particular particles (correctly handels the
  // case of particles with non-zero width)
  double threshold2(const particle& vm, const particle& recoil) const;

  // utility function
  lp_gamma_event make_event(const lp_gamma_data& initial, const double t,
                            const double b, particle vm1, particle X1,
                            const double xs, const double R);

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;

  // scattering amplitude
  class oleksii_2vmp_amplitude* ampl_;
  class oleksii_2vmp_slope* slope_;

  // cross section settings
  const double T0_;       // subtraction constant for dispersion relation
  const double R_vm_c_;   // c-parameter for R
  const double R_vm_n_;   // n-parameter for R
  const double dipole_n_; // n-parameter for dipole factor

  // t-range and cross setion maxima
  const double max_b_; // upper limit to b parameter
  const interval<double> max_t_range_;
  const double max_;
};

} // namespace lp_gamma
} // namespace pcsim

#endif
