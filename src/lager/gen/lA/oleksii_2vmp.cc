// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2020 Sylvester Joosten <sjoosten@anl.gov>
// 
// This file is part of lAger.
// 
// lAger is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Shoftware Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// lAger is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with lAger.  If not, see <https://www.gnu.org/licenses/>.
// 

#include "oleksii_2vmp.hh"
#include <TMath.h>
#include <lager/core/logger.hh>
#include <lager/physics/kinematics.hh>
#include <lager/physics/photon.hh>
#include <lager/physics/vm.hh>

#include <Math/RootFinderAlgorithms.h>
#include <Math/WrappedTF1.h>
#include <TF1.h>
#include <cmath>
#include <complex>
#include <gsl/gsl_integration.h>

namespace lager {
namespace lA {
// local constants etc.
namespace {
// unit conversion
// Planks constant (eV * S, CODATA)
const double khbar = 4.135667662e-15 / TMath::TwoPi();
// Speed of light (m/s)
const double kc = 299792458.;
// hhbarc (fm * GeV)
const double khbarc = khbar * kc * 1e6;
// khbarc2 (nb * GeV^2, using 1fm^2 = 10mb^2 =1e7nb^2)
const double khbarc2 = khbarc * khbarc * 1e7;

// VM decay constant (in GeV)
const double kfu = 0.238;
const double kfj = 0.278;
// electric charge
const double kalpha = 1 / 137.;
const double ke = sqrt(4 * TMath::Pi() * kalpha);
// Masses (in GeV)
const double kMu = 9.46030;
const double kMj = 3.096916;
const double kMp = 0.938272;
const double kMu2 = kMu * kMu;
const double kMp2 = kMp * kMp;
const double kMj2 = kMj * kMj;
// photon 3-momentum in CM frame
auto qgp = [](const double s) { return (s - kMp2) / (2. * sqrt(s)); };
// Vector meson 3-momentum in CM frame
auto qvp = [](const double s, const double Mv) {
  const double Mv2 = Mv * Mv;
  return sqrt(1 / (4 * s) * (s - (Mv2 + 2 * Mv * kMp + kMp2)) *
              (s - (Mv2 - 2 * Mv * kMp + kMp2)));
};
// fit constants
const double kC_el_jpsi = .1;
const double kC_in_jpsi = 20.51;
const double kC_el_upsilon = 14.30e-3;
const double kC_in_upsilon = 18.95;
const double kb_el = 1.27;
const double kb_in = 3.53;
const double ka_el = 1.38;
const double ka_in = 1.2;
// crossing variable v
auto v = [](const double s, const double Mv) {
  return 0.5 * (s - kMp2 - Mv * Mv);
};
// crossing variables at threshold
const double kv_el_jpsi = kMp * kMj;
const double kv_el_upsilon = kMp * kMu;
const double kv_in_jpsi = 5.66;        // at DD threshold
const double kv_in_upsilon = 20.89999; // at BB threshold
} // namespace
// amplitude and cross section implementation
class oleksii_2vmp_amplitude {
  struct IntegrandParam {
    const double nu;
    const oleksii_2vmp_amplitude* that;
  };

public:
  oleksii_2vmp_amplitude(pdg_id id, double subtraction_constant)
      : T0_{subtraction_constant}
      , v_el_{id == pdg_id::upsilon ? kv_el_upsilon : kv_el_jpsi}
      , v_in_{id == pdg_id::upsilon ? kv_in_upsilon : kv_el_jpsi}
      , C_el_{id == pdg_id::upsilon ? kC_el_upsilon : kC_el_jpsi}
      , C_in_{id == pdg_id::upsilon ? kC_in_upsilon : kC_in_jpsi}
      , Mv_{id == pdg_id::upsilon ? kMu : kMj}
      , fv_{id == pdg_id::upsilon ? kfu : kfj}
      , w_{gsl_integration_workspace_alloc(2048)} {
    F_.function = [](double nup, void* param) {
      if (!param) {
        return 0.;
      }
      const IntegrandParam* p = static_cast<IntegrandParam*>(param);
      const double nu = p->nu;
      const oleksii_2vmp_amplitude* that = p->that;
      const double result = (that->ImT_nu(nup) / nup - that->ImT_nu(nu) / nu) /
                            (nup * nup - nu * nu);
      return result;
    };
  }
  ~oleksii_2vmp_amplitude() { gsl_integration_workspace_free(w_); }
  // Im(T) as a function of nu
  double ImT_nu(const double nu) const {
    double ImT = 0;
    if (nu > v_el_) {
      ImT += C_el_ * pow(1 - v_el_ / nu, kb_el) * pow(nu / v_el_, ka_el);
    }
    if (nu > v_in_) {
      ImT += C_in_ * pow(1 - v_in_ / nu, kb_in) * pow(nu / v_in_, ka_in);
    }
    return ImT;
  }
  double ImT(const double W) const {
    const double nu = v(W * W, Mv_);
    return ImT_nu(nu);
  }
  double ReT(const double W) const {
    const double nu = v(W * W, Mv_);
    return T0_ + 2 * TMath::InvPi() * nu * nu * DispersionIntegral(nu);
  }
  double T0() const { return T0_; }
  double Mv() const { return Mv_; }
  std::complex<double> T(const double W) const { return {ReT(W), ImT(W)}; }
  // differential cross section at t=0 (nb/GeV^2)
  double dsdt_t0(const double W) const {
    const double s = W * W;
    const double T2 = std::norm(T(W));
    const double factor =
        std::pow(sqrt(4. * TMath::Pi() * kalpha) * fv_ / Mv_, 2) * 1 /
        (64 * TMath::Pi() * s * qgp(s) * qgp(s));
    return factor * T2 * khbarc2;
  }
  // total integrated photo-production cross section (nb)
  double sigma_el(const double W) const {
    const double s = W * W;
    const double rt_s = W;
    if (W < (kMp + Mv_)) {
      return 0.;
    }
    return std::pow(ke * fv_ / Mv_, 2) * 1 / (2 * rt_s * qgp(s)) *
           (qvp(s, Mv_) / qgp(s)) * C_el_ *
           std::pow(1 - v_el_ / v(s, Mv_), kb_el) *
           std::pow(v(s, Mv_) / v_el_, ka_el) * khbarc2;
  }

private:
  double DispersionIntegral(const double nu) const {
    if (nu < v_el_) {
      return 0;
    }
    double result, error;
    IntegrandParam param{nu, this};
    F_.params = &param;
    gsl_integration_qagiu(&F_, v_el_, 0, 1e-5, 2048, w_, &result, &error);
    return result +
           ImT_nu(nu) / nu * std::log(std::fabs((v_el_ + nu) / (v_el_ - nu))) /
               (2 * nu);
  }
  const double T0_;
  const double v_el_;
  const double v_in_;
  const double C_el_;
  const double C_in_;
  const double Mv_;
  const double fv_;
  mutable gsl_integration_workspace* w_;
  mutable gsl_function F_;
};
// B slope implementation
class oleksii_2vmp_slope {
public:
  oleksii_2vmp_slope(const oleksii_2vmp_amplitude& ampl)
      : ampl_{ampl}
      , equation_{"B_equation",
                  [&](double* bb, double* par) {
                    const double B = bb[0];
                    const double W = par[0];
                    const double Q2 = par[1];
                    const auto tlim =
                        physics::t_range(W * W, Q2, kMp, ampl.Mv(), kMp);
                    return B -
                           this->B0(W) *
                               (exp(B * tlim.max) - exp(B * tlim.min));
                  },
                  0, 10, 2}
      , fwrap_{equation_} {}

  double B0(const double W) const {
    const double sigma = ampl_.sigma_el(W);
    if (sigma > 0) {
      const double result = ampl_.dsdt_t0(W) / sigma;
      return result;
    }
    return 0;
  }
  double B(const double W, const double Q2) const {
    equation_.SetParameters(W, Q2);
    brf_.SetFunction(fwrap_, 1e-7, 10);
    brf_.Solve();
    const double result = brf_.Root();
    return result;
  }

private:
  const oleksii_2vmp_amplitude& ampl_;
  mutable ROOT::Math::Roots::Brent brf_;
  mutable TF1 equation_;
  mutable ROOT::Math::WrappedTF1 fwrap_;
};

// =============================================================================
// Constructor for lA::oleksii_2vmp
// =============================================================================
oleksii_2vmp::oleksii_2vmp(const configuration& cf, const string_path& path,
                           std::shared_ptr<TRandom> r)
    : base_type{r}
    , recoil_{pdg_id::p}
    , vm_{cf.get<std::string>(path / "vm_type")}
    , ampl_{new oleksii_2vmp_amplitude(vm_.type(), cf.get<double>(path / "T0"))}
    , slope_{new oleksii_2vmp_slope(*ampl_)}
    , T0_{cf.get<double>(path / "T0")}
    , R_vm_c_{cf.get<double>(path / "R_vm_c")}
    , R_vm_n_{cf.get<double>(path / "R_vm_n")}
    , dipole_n_{cf.get<double>(path / "dipole_n")}
    , max_b_range_{calc_max_b_range(cf)}
    , max_t_range_{calc_max_t_range(cf)}
    , max_exp_b0t_range_{exp(max_b_range_.max * max_t_range_.min), 1.}
    , max_{calc_max_xsec(cf)} {
  LOG_INFO("oleksii_2vmp",
           "t range [GeV^2]: [" + std::to_string(max_t_range_.min) + ", " +
               std::to_string(max_t_range_.max) + "]");
  LOG_INFO("oleksii_2vmp",
           "b parameter lower limit [1/GeV^2]: " +
               std::to_string(max_b_range_.min));
  LOG_INFO("oleksii_2vmp",
           "b parameter upper limit [1/GeV^2]: " +
               std::to_string(max_b_range_.max));
  LOG_INFO("oleksii_2vmp", "subtraction constant T0: " + std::to_string(T0_));
  LOG_INFO("oleksii_2vmp", "R_vm c-parameter: " + std::to_string(R_vm_c_));
  LOG_INFO("oleksii_2vmp",
           "R_vm n-parameter (power): " + std::to_string(R_vm_n_));
  LOG_INFO("oleksii_2vmp", "'Dipole' FF power: " + std::to_string(dipole_n_));
  LOG_INFO("oleksii_2vmp", "VM: " + std::string(vm_.pdg()->GetName()));
  LOG_INFO("oleksii_2vmp", "recoil: " + std::string(recoil_.pdg()->GetName()));
}

lA_event oleksii_2vmp::generate(const lA_data& initial) {

  // generate a mass() in case of non-zero width, initialize the particles
  particle vm = {vm_.type(), rng()};
  particle recoil = {recoil_.type(), rng()};

  // shortcuts
  const auto& gamma = initial.photon();
  const auto& target = initial.target();

  // check if enough energy available
  if (gamma.W2() < threshold2(vm, recoil)) {
    LOG_JUNK(
        "oleksii_2vmp",
        "Not enough phase space available - W2: " + std::to_string(gamma.W2()) +
            " < " + std::to_string(threshold2(vm, recoil)));
    return lA_event{0.};
  }

  // generate a phase space point
  const double t =
      std::log(rng()->Uniform(max_exp_b0t_range_.min, max_exp_b0t_range_.max)) /
      max_b_range_.min;
  const double b = slope_->B(gamma.W(), gamma.Q2());

  LOG_JUNK("oleksii_2vmp",
           "t: " + std::to_string(t) + ", b: " + std::to_string(b));

  // check if kinematically allowed
  if (physics::t_range(gamma.W2(), gamma.Q2(), target.particle().mass(),
                       vm.mass(), recoil.mass())
          .excludes(t)) {
    LOG_JUNK("oleksii_2vmp", "t outside of the allowed range for this W2")
    return lA_event{0.};
  }

  // evaluate the cross section
  const double xs_R = R(gamma.Q2());
  const double xs_dipole = dipole(gamma.Q2());
  const double xs_photo = dsigma_dexp_b0t(gamma.W2(), t, b);
  const double xs = (1 + gamma.epsilon() * xs_R) * xs_dipole * xs_photo;

  LOG_JUNK("oleksii_2vmp",
           "xsec: " + std::to_string(xs_photo) + " < " + std::to_string(max_));
  LOG_JUNK("oleksii_2vmp", "R: " + std::to_string(xs_R));
  LOG_JUNK("oleksii_2vmp", "dipole: " + std::to_string(xs_dipole));

  // return a new VM event
  return make_event(initial, t, b, vm, recoil, xs, xs_R);
}

interval<double> oleksii_2vmp::calc_max_b_range(const configuration& cf) const {
  return {calc_min_b(), calc_max_b(cf)};
}
double oleksii_2vmp::calc_max_b(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the
  // lepton beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                        cf.get<double>("beam/lepton/energy")};
  const particle target{cf.get<std::string>("beam/ion/particle_type"),
                        cf.get_vector3<particle::XYZVector>("beam/ion/dir"),
                        cf.get<double>("beam/ion/energy")};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W
  const double Wmax =
      opt_W_range ? fmin(opt_W_range->max, (photon.p() + target.p()).M())
                  : (photon.p() + target.p()).M();
  // get our b
  const double result = slope_->B(Wmax, 0);
  std::cout << result << std::endl;
  return result;
}
// Get smallest possible b value (in case of no binding && at threshold)
double oleksii_2vmp::calc_min_b() const {
  oleksii_2vmp_amplitude nb_ampl{vm_.type(), 0};
  oleksii_2vmp_slope nb_slope{nb_ampl};
  return nb_slope.B((vm_.mass() + kMp) * 1.05, 0);
}

// =============================================================================
// oleksii_2vmp::calc_max_xsec(cf)
//
// Utility function for the generator initialization
//
// max cross section as defined by the beam/target and phase-space settings
//
// =============================================================================
double oleksii_2vmp::calc_max_xsec(const configuration& cf) const {
  // get the extreme beam/lepton parameters (where the photon carries all of the
  // lepton beam/lepton energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                        cf.get<double>("beam/lepton/energy")};
  const particle target{cf.get<std::string>("beam/ion/particle_type"),
                        cf.get_vector3<particle::XYZVector>("beam/ion/dir"),
                        cf.get<double>("beam/ion/energy")};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W
  const double W2max = opt_W_range ? fmin(opt_W_range->max * opt_W_range->max,
                                          (photon.p() + target.p()).M2())
                                   : (photon.p() + target.p()).M2();
  // possible max at t=tmin and max W (where b is maximum)
  const double max = dsigma_dexp_b0t(W2max, max_t_range_.max, max_b_range_.min);
  // other possible maximum near minimum b
  // return max * 4; // need pretty generous padding to not exceed this
  // "maximum"
  //                // (there are higher maxima near lower W where b is smaller)
  return max;
}

// =============================================================================
// oleksii_2vmp::calc_max_t_range(cf)
//
// Utility function for the generator initialization
//
// max t-range occurs for:
//  * photon that carries all of the beam energy (or when we have reached the
//    user-defined maximum value of W max)
//  * maximum Q2 for the given W (or zero for real photons) for the lower t
//  bound
//  * Q2 = 0 for the upper t bound (tmin)
// =============================================================================
interval<double> oleksii_2vmp::calc_max_t_range(const configuration& cf) const {
  // get the extreme beam parameters (where the photon carries all of the
  // lepton beam energy
  const particle photon{pdg_id::gamma,
                        cf.get_vector3<particle::XYZVector>("beam/lepton/dir"),
                        cf.get<double>("beam/lepton/energy")};
  const particle target{cf.get<std::string>("beam/ion/particle_type"),
                        cf.get_vector3<particle::XYZVector>("beam/ion/dir"),
                        cf.get<double>("beam/ion/energy")};
  // check if we have a user-defined W-range set
  const auto opt_W_range = cf.get_optional_range<double>("photon/W_range");
  // get the maximum W (where the t-range is the largest)
  const double W2max = opt_W_range ? fmin(opt_W_range->max * opt_W_range->max,
                                          (photon.p() + target.p()).M2())
                                   : (photon.p() + target.p()).M2();
  // some shortcuts
  const double M_nu = (photon.p()).Dot(target.p());
  const double Q2max = target.mass2() + 2 * M_nu - W2max;
  const double Q2min = 0.;
  // calculate the corresponding t range
  // in case of particles with non-zero width, we use M - 4 x sigma
  //
  // tlim1: the correct lower t limit (tmax)
  const auto tlim1 =
      physics::t_range(W2max, (Q2max > 1e-10 ? Q2max : 0.), target.mass(),
                       vm_.pole_mass() - vm_.width() * 4.,
                       recoil_.pole_mass() - recoil_.width() * 4);
  // tlim2: the correct upper t limit (tmin)
  const auto tlim2 = physics::t_range(
      W2max, 0, target.mass(), vm_.pole_mass() - vm_.width() * 4.,
      recoil_.pole_mass() - recoil_.width() * 4);
  const auto tlim = interval<double>(tlim1.min, tlim2.max);
  return tlim;
}
// =============================================================================
// oleksii_2vmp::dsigma_dt()
// oleksii_2vmp::jacobian()
// oleksii_2vmp::R()
// oleksii_2vmp::dipole()
//
// Utility functions to calculate the cross section components
// =============================================================================
double oleksii_2vmp::dsigma_dexp_b0t(const double W2, const double t,
                                     const double b) const {
  // xsec at t=0
  const double t0 = ampl_->dsdt_t0(sqrt(W2));
  // form factor
  const double ff = exp(b * t);
  // jacobian for dt -> d(b0t)
  double jacobian = 1 / max_b_range_.min;
  // jacobian for d(b0t) -> d(exp_b0t)
  jacobian *= exp(-max_b_range_.min * t);
  return t0 * ff * jacobian;
}

double oleksii_2vmp::jacobian(const double t) const {
  return max_b_range_.min * exp(max_b_range_.min * t);
}
double oleksii_2vmp::R(const double Q2) const {
  return physics::R_vm_martynov(Q2, vm_.mass(), R_vm_c_, R_vm_n_);
}
double oleksii_2vmp::dipole(const double Q2) const {
  return physics::dipole_ff_vm(Q2, vm_.mass(), dipole_n_);
}
// =============================================================================
// oleksii_2vmp::threshold2()
//
// utility function returns the production threshold squared for the chosen
// particles. Important to re-calculate in case of a particle with non-zero
// widht.
// =============================================================================
double oleksii_2vmp::threshold2(const particle& vm,
                                const particle& recoil) const {

  return recoil.mass2() + vm.mass2() + 2 * vm.mass() * recoil.mass();
}

// =============================================================================
// create the lA_event dataf, calculates the final state four-vectors in
// the lab-frame
// =============================================================================
lA_event oleksii_2vmp::make_event(const lA_data& initial, const double t,
                                  const double b, particle vm, particle X,
                                  const double xs, const double R) {
  const auto& gamma = initial.photon();
  const auto& target = initial.target();

  lA_event e{initial, xs, 1., R};

  e.update_jacobian(jacobian(t));

  // utility shortcuts
  const double W2 = gamma.W2();
  const double W = sqrt(W2);
  const double Q2 = gamma.Q2();
  const double y = gamma.y();
  const double nu = gamma.nu();
  const double Mt2 = target.particle().mass2();
  const double Mr2 = X.mass2();
  const double Mv2 = vm.mass2();

  // create our final state particles in the CM frame
  // energies and momenta
  const double Et_cm = (W2 + Q2 + Mt2) / (2. * W);
  const double Pt_cm = sqrt(Et_cm * Et_cm - Mt2);
  const double Er_cm = (W2 - Mv2 + Mr2) / (2. * W);
  const double Pr_cm = sqrt(Er_cm * Er_cm - Mr2);
  const double Ev_cm = (W2 + Mv2 - Mr2) / (2. * W);
  const double Pv_cm = sqrt(Ev_cm * Ev_cm - Mv2);

  // get the VM and recoil CM 4-vectors
  // NOTE:
  //  * theta is the change in angle from the initial state, don't forget the
  //    target originally was flying backwards (already has theta=pi)
  // calculate the scattering angle theta
  const double ctheta_cm =
      (t + 2 * Et_cm * Er_cm - Mt2 - Mr2) / (2 * Pt_cm * Pr_cm);
  const double theta_cm = std::acos(ctheta_cm);
  const double phi_cm = rng()->Uniform(0, TMath::TwoPi());

  // create the CM 4-vectors
  const particle::Polar3DVector p3_v{Pv_cm, theta_cm, phi_cm};
  vm.p() = {p3_v.X(), p3_v.Y(), p3_v.Z(), Ev_cm};
  const particle::Polar3DVector p3_r{Pr_cm, theta_cm + TMath::Pi(), phi_cm};
  X.p() = {p3_r.X(), p3_r.Y(), p3_r.Z(), Er_cm};

  // calculate some necessary boost and rotation vectors
  particle targ = target.particle();
  particle phot = gamma.particle();
  // lab frame to target rest frame (trf)
  particle::Boost boost_to_trf{targ.p().BoostToCM()};
  targ.boost(boost_to_trf);
  phot.boost(boost_to_trf);
  // describe photon in a target rest frame where the photon moves along the
  // z-axis
  particle::XYZTVector p_phot_z{0, 0, phot.momentum(), phot.energy()};
  // boost to CM frame from this rotated trf
  particle::Boost boost_from_cm{-((targ.p() + p_phot_z).BoostToCM())};

  // go to lab frame
  // 1. CM -> rotated TRF
  vm.boost(boost_from_cm);
  X.boost(boost_from_cm);
  // 2. rotated TRF -> TRF
  vm.rotate_uz(phot.p());
  X.rotate_uz(phot.p());
  // 3. TRF -> lab
  vm.boost(boost_to_trf.Inverse());
  X.boost(boost_to_trf.Inverse());

  // update VM particle status
  vm.update_status(particle::status_code::UNSTABLE_SCHC);

  // set vertex info
  vm.vertex() = phot.vertex();
  X.vertex() = phot.vertex();

  // add to the event
  e.add_leading(vm, e.photon_index(), e.target_index());
  e.add_recoil(X, e.photon_index(), e.target_index());

  // all done!
  return e;
}

} // namespace lA
} // namespace lager
