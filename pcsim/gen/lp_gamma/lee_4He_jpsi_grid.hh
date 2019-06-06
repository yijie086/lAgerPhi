#ifndef PCSIM_GEN_LP_GAMMA_LEE_4HE_JPSI_GRID_LOADED
#define PCSIM_GEN_LP_GAMMA_LEE_4HE_JPSI_GRID_LOADED

#include <TGraph2D.h>
#include <memory>
#include <pcsim/core/generator.hh>
#include <pcsim/core/particle.hh>
#include <pcsim/gen/lp_gamma/generator.hh>
#include <pcsim/gen/lp_gamma_event.hh>

namespace pcsim {
namespace lp_gamma {

// =============================================================================
// lp_gamma::lee_4He_jpsi_grid
//
// gamma + p -> VM + X process
//
// Uses the following expressions (cf. pcsim/physics/vm.hh)
//  * R (sigma_L/sigma_T):
//        R_vm_martynov(...)
//  * Dipole FF for sigma_gamma -> sigma_t:
//        dipole_ff_vm_hermes(...)
//  * t-channel cross section:
//        dsigma_dexp_bt_brodsky(...)
// =============================================================================
class lee_4He_jpsi_grid : public lp_gamma::generator {
public:
  using base_type = lp_gamma::generator;

  lee_4He_jpsi_grid(const configuration& cf, const string_path& path,
                    std::shared_ptr<TRandom> r);
  virtual lp_gamma_event generate(const lp_gamma_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return max_t_range_.width(); }

private:
  double calc_max_xsec(const configuration& cf) const;
  interval<double> calc_max_t_range(const configuration& cf) const;

  // cross section component evaluation
  double dsigma_dt(const double W2, const double t, const double Mt) const;
  double R(const double Q2) const;
  double dipole(const double Q2) const;

  // jacobian for d/dexp_bt -> d/dt
  double jacobian(const double t) const;

  // threshold squared for these particular particles (correctly handels the
  // case of particles with non-zero width)
  double threshold2(const particle& vm, const particle& recoil) const;

  // utility function
  lp_gamma_event make_event(const lp_gamma_data& initial, const double t,
                            particle vm1, particle X1, const double xs,
                            const double R);

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;

  // cross section settings
  std::unique_ptr<TGraph2D> grid_;
  const double R_vm_c_;   // c-parameter for R
  const double R_vm_n_;   // n-parameter for R
  const double dipole_n_; // n-parameter for dipole factor

  // t-range and cross setion maxima
  const interval<double> max_t_range_;
  const double max_;
};

} // namespace lp_gamma
} // namespace pcsim

#endif
