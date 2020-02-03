#ifndef LIEGE_GEN_LP_GAMMA_OLEKSII_JPSI_BH_LOADED
#define LIEGE_GEN_LP_GAMMA_OLEKSII_JPSI_BH_LOADED

#include <liege/core/generator.hh>
#include <liege/core/particle.hh>
#include <liege/gen/lp_gamma/generator.hh>
#include <liege/gen/lp_gamma_event.hh>

namespace liege {
namespace lp_gamma {

// =============================================================================
// lp_gamma::oleksii_jpsi_bh
//
// gamma + p -> VM + X process
//
// Uses the following expressions (cf. liege/physics/vm.hh)
//  * R (sigma_L/sigma_T):
//        R_vm_martynov(...)
//  * Dipole FF for sigma_gamma -> sigma_t:
//        dipole_ff_vm_hermes(...)
//  * t-channel cross section:
//        dsigma_dexp_bt_brodsky(...)
// =============================================================================
class oleksii_jpsi_bh : public lp_gamma::generator {
public:
  using base_type = lp_gamma::generator;

  oleksii_jpsi_bh(const configuration& cf, const string_path& path,
                  std::shared_ptr<TRandom> r);
  virtual lp_gamma_event generate(const lp_gamma_data&);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const {
    // last factor is the cos(theta) range
    // return max_t_range_.width() *
    return (exp(1.13 * max_t_range_.max) - exp(1.13 * max_t_range_.min)) *
           (Mll_range_.max * Mll_range_.max - Mll_range_.min * Mll_range_.min) *
           TMath::TwoPi() * 2.;
  }

private:
  interval<double> calc_max_t_range(const configuration& cf) const;

  // threshold squared for these particular particles (correctly handels the
  // case of particles with non-zero width)
  double threshold2(const particle& vm, const particle& recoil) const;

  // utility function
  lp_gamma_event make_event(const lp_gamma_data& initial, const double t,
                            particle vm1, particle X1, const double xs,
                            const double thetaCM_el,
                            const double phiCM_el) const;

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;

  // t-range and cross setion maxima
  const interval<double> max_t_range_;
  const interval<double> Mll_range_;
  const double max_;

  // subtraction constant
  const double T_0_;

  // theta acceptance to cut out colinear enhancements
  const interval<double> p_range_;
  const interval<double> theta_range_;
};

} // namespace lp_gamma
} // namespace liege

#endif
