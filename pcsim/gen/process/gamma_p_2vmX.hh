#ifndef PCSIM_GEN_PROCESS_GAMMA_P_2VMX_LOADED
#define PCSIM_GEN_PROCESS_GAMMA_P_2VMX_LOADED
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
// gamma_p_2vmX_data
//
// VM production data
// =============================================================================
class gamma_p_2vmX_data : public generator_data {
public:
  gamma_p_2vmX_data() = default;
  gamma_p_2vmX_data(const gamma_p_2vmX_data&) = default;
  gamma_p_2vmX_data& operator=(const gamma_p_2vmX_data&) = default;

  gamma_p_2vmX_data(const double xs) : generator_data{xs} {}

  gamma_p_2vmX_data(const beam::photon_data& photon, const beam::data& target,
                    const double t, const particle& vm1, const particle& X1,
                    const double xs, const double phi, const double R = 0);

  double t() const { return t_; }
  double xv() const { return xv_; }
  double Q2plusMv2() const { return Q2plusMv2_; }
  double R() const { return R_; }
  const particle& vm() { return vm_; }
  const particle& X() { return X_; };

private:
  double t_;
  double xv_;
  double Q2plusMv2_;
  double R_;
  particle vm_;
  particle X_;
};

// =============================================================================
// process::gamma_p_2vmX
//
// abstract base class for all gamma_p_2vmX processes
// =============================================================================
class gamma_p_2vmX
    : public generator<gamma_p_2vmX_data, beam::photon_data, beam::data> {
public:
  static factory<gamma_p_2vmX> factory;

  gamma_p_2vmX(std::shared_ptr<TRandom> r) : generator{std::move(r)} {}
};

// =============================================================================
// process::gamma_p_2vmX_brodksy
//
// gamma + p -> VM + X process
//
// Uses the following expressions (cf. pcsim/physics/vm.hh)
//  * R (sigma_L/sigma_T):
//        R_vm_martynov(...)
//  * Dipole FF for sigma_gamma -> sigma_t:
//        dipole_ff_vm_hermes(...)
//  * t-channel cross section:
//        dsigma_dexp_bt_brodksy(...)
// =============================================================================
class gamma_p_2vmX_brodksy : public gamma_p_2vmX {
public:
  gamma_p_2vmX_brodsky(const configuration& cf, const string_path& path,
                       std::shared_ptr<TRandom> r);
  virtual gamma_p_2vmX_data generate(const beam::photon_data& photon,
                                     const beam::data& target);
  virtual double max_cross_section() const { return max_; }
  virtual double phase_space() const { return max_exp_bt_range_.width(); }

private:
  double calc_max_xsec(const configuration& cf) const;
  interval<double> calc_max_t_range(const configuration& cf) const;

  // kinematics range
  interval<double> exp_bt_range(const double W2, const double Q2, const double Mt) const;

  // cross section component evaluation
  double dsigma_dexp_bt(const double W2, const double Mt) const;
  double R(const double Q2) const;
  double dipole(const double Q2) const;

  // recoil and vm particle info
  const particle recoil_;
  const particle vm_;
  const double threshold2_; // threshold squared

  // cross section settings
  const double photo_b_;   // target FF constant
  const double photo_c2g_; // 2-gluon amplitude norm
  const double photo_c3g_; // 3-gluon amplitude norm
  const double R_vm_c_;    // c-parameter for R
  const double R_vm_n_;    // n-parameter for R
  const double dipole_n_;  // n-parameter for dipole factor

  // t-range and cross setion maxima
  const interval<double> max_t_range_;
  const interval<double> max_exp_bt_range_;
  const double max_;
};

} // namespace process
} // namespace pcsim

#endif
