// ROOT module with TF1 versions of various physics functions
// NOTE: requires that setup_root.py was loaded prior to compilation,
//       or for the relevant modules to already be manually loaded

#include "physics.hh"
#include <Math/Vector4D.h>
#include <TF1.h>
#include <algorithm>
#include <cmath>
#include <pcsim/physics/kinematics.hh>
#include <pcsim/physics/photon.hh>
#include <pcsim/physics/vm.hh>

namespace vm {
// epsilon as a function of Q2 (W and and beam energy are parameters)
// x[0]: Q2
// par[1]: W
// par[2]: Ebeam
TF1* epsilon_Q2(double Q2min, double Q2max) {
  return new TF1("epsilon_Q2_WE",
                 [](double* xx, double* par) {
                   const double Q2 = xx[0];
                   const double W = par[0];
                   const double Ebeam = par[1];
                   return epsilon(xx[0], par[0], par[1]);
                 },
                 Q2min, Q2max, 2);
}

// =============================================================================
// t-channel photo-production cross section for heavy vector mesons
// Values are in units of nb/GeV^2.
//
// Formulism from
//      brodsky et. al., Phys.Lett.B498:23-28,2001
//      (http://arxiv.org/abs/hep-ph/0010343)
//
// =============================================================================
// x[0]: t (sign force-added, in case we want to draw -t)
// par[0]: 2-gluon amplitude
// par[1]: 3-gluon amplitude
// par[2]: b
// par[3]: VM mass
// par[4]: s/W^2
TF1* dsigma_dt_vm_brodsky(const double tmin, const double tmax) {
  return new TF1("dsigma_dt_vm_brodsky",
                 [](double* xx, double* par) {
                   const double Mt = M_P;
                   const double t = -std::fabs(xx[0]);
                   const double s = par[4];
                   const double Mv = par[3];
                   const double b = par[2];
                   const double c3g = par[1];
                   const double c2g = par[0];
                   return pcsim::physics::dsigma_dt_vm_brodsky(s, t, Mt, Mv, b,
                                                               c2g, c3g);
                 },
                 std::min(tmin, tmax), std::max(tmin, tmax), 5);
}
// same but differential in d/d(exp(bt)) for more efficient integration
// x[0]: t (unused)
// par[0]: 2-gluon amplitude
// par[1]: 3-gluon amplitude
// par[2]: b
// par[3]: VM mass
// par[4]: s/W^2
TF1* dsigma_dexp_bt_vm_brodsky(double tmin, double tmax) {
  return new TF1("dsigma_dexp_bt_vm_brodsky",
                 [](double* xx, double* par) {
                   const double Mt = M_P;
                   const double t = -std::fabs(xx[0]);
                   const double s = par[4];
                   const double Mv = par[3];
                   const double b = par[2];
                   const double c3g = par[1];
                   const double c2g = par[0];
                   return pcsim::physics::dsigma_dexp_bt_vm_brodsky(
                       s, Mt, Mv, b, c2g, c3g);
                 },
                 std::min(tmin, tmax), std::max(tmin, tmax), 5);
}
// t-integrated cross section as a function of W
// x[0]: W
// par[0]: 2-gluon Amplitude
// par[1]: 3-gluon Amplitude
// par[2]: b
// par[3]: VM mass
TF1* sigma_vm_brodsky_W(double Wmin, double Wmax) {
  return new TF1("sigma_vm_brodsky_W",
                 [](double* xx, double* par) {
                   const double s = xx[0] * xx[0];
                   const double Mv = par[3];
                   const double b = par[2];
                   const double c3g = par[1];
                   const double c2g = par[0];
                   // get the t-range, transform to exp(bt)
                   auto tl = pcsim::physics::t_range(s, 0., M_P, Mv, M_P);
                   tl.min = std::exp(b * tl.min);
                   tl.max = std::exp(b * tl.max);
                   // We integrate the dsigma/dexp(bt) formula, which is
                   // constant as a function of t. Hence, we only need to
                   // evaluate the function once, and multiply with the
                   // exp(bt)-range
                   const double delta_t = tl.min - tl.max;
                   double xs[1] = {tl.min};
                   const double integral =
                       pcsim::physics::dsigma_dexp_bt_vm_brodsky(s, M_P, Mv, b,
                                                                 c2g, c3g) *
                       delta_t;
                   return integral;
                 },
                 Wmin, Wmax, 4.);
}
// =============================================================================
// R VM parameterization
//
// from
//      Martynov, et. al., “Photoproduction of Vector Mesons in the Soft Dipole
//      Pomeron Model.” PRD 67 (7), 2003. doi:10.1103/PhysRevD.67.074023.
// eq. 31.
//
// J/psi parameters c(or a) and n
//  * c: 2.164
//  * n: 2.131
// from eq. 18:
//
//      R. Fiore et al., "Exclusive Jpsi electroproduction in a dual model."
//      PRD80:116001, 2009"
// =============================================================================
// x[0]: Q2
// par[0]: Mv
// par[1]: c
// par[2]: n
TF1* R_vm_martynov(const double Q2min, const double Q2max) {
  return new TF1("R_vm_martynov",
                 [](double* xx, double* par) {
                   const double Q2 = xx[0];
                   const double Mv = par[0];
                   const double c = par[1];
                   const double n = par[2];
                   return pcsim::physics::R_vm_martynov(Q2, Mv, c, n);
                 },
                 Q2min, Q2max, 3);
}
// =============================================================================
// Dipole form factor to relate real photo-production cross section to sigma_T
//
// This form deviates from the classical VMD form (which has a fixed power of
// 2), to better fit the world data for rho0 production.
//
// Cf.:
//      A. Airapetian et al, "Exclusive Leptoproduction of rho0 Mesons on
//      Hydrogen at Intermediate W Values", EPJ C 17 (2000) 389-398
//
//      Adams et al., "Diffractive production of ρ0 mesons in muon–proton
//      interactions 470 GeV", ZPC74 (1997) 237-261.
// =============================================================================
// x[0]: Q2
// par[0]: Mv
// par[1]: n
TF1* dipole_ff_vm(const double Q2min, const double Q2max) {
  return new TF1("dipole_ff_vm",
                 [](double* xx, double* par) {
                   const double Q2 = xx[0];
                   const double Mv = par[0];
                   const double n = par[1];
                   return pcsim::physics::dipole_ff_vm(Q2, Mv, n);
                 },
                 Q2min, Q2max, 2);
}

// =============================================================================
// t-integrated VM cross section as a function of (Q2, W) (x-val: W)
// =============================================================================
// x[0]: W
// par[0]: 2-gluon Amplitude
// par[1]: 3-gluon Amplitude
// par[2]: b
// par[3]: VM mass
// par[4]: Q2
// par[5]: W (unused)
// par[6]: Ebeam
// par[7]: dipole_n
// par[8]: R_c
// par[9]: R_n
TF1* sigma_vm_brodsky_W_Q2(const double Wmin, const double Wmax) {
  return new TF1(
      "sigma_vm_brodsky_W_Q2",
      [](double* xx, double* par) {
        const double s = xx[0] * xx[0];
        const double W = xx[0];
        const double R_n = par[9];
        const double R_c = par[8];
        const double dipole_n = par[7];
        const double Ebeam = par[6];
        // par[5] unused
        const double Q2 = par[4];
        const double Mv = par[3];
        const double b = par[2];
        const double c3g = par[1];
        const double c2g = par[0];
        // get the t-range, transform to exp(bt)
        auto tl = pcsim::physics::t_range(s, Q2, M_P, Mv, M_P);
        tl.min = std::exp(b * tl.min);
        tl.max = std::exp(b * tl.max);
        // We integrate the dsigma/dexp(bt) formula, which is
        // constant as a function of t. Hence, we only need to
        // evaluate the function once, and multiply with the
        // exp(bt)-range
        const double delta_t = tl.min - tl.max;
        double xs[1] = {tl.min};
        const double integral =
            pcsim::physics::dsigma_dexp_bt_vm_brodsky(s, M_P, Mv, b, c2g, c3g) *
            delta_t;
        // add rest of Q2 dependence
        return integral * pcsim::physics::dipole_ff_vm(Q2, Mv, dipole_n) *
               (1. + pcsim::physics::R_vm_martynov(Q2, Mv, R_c, R_n) *
                         epsilon(Q2, W, Ebeam));
      },
      Wmin, Wmax, 10);
}

// =============================================================================
// VM cross section differential in t as a function of (Q2, W) (x-val: t)
// =============================================================================
// x[0]: t
// par[0]: 2-gluon Amplitude
// par[1]: 3-gluon Amplitude
// par[2]: b
// par[3]: VM mass
// par[4]: Q2
// par[5]: W
// par[6]: Ebeam
// par[7]: dipole_n
// par[8]: R_c
// par[9]: R_n
TF1* dsigma_dt_vm_brodsky_Q2W(const double tmin, const double tmax) {
  return new TF1(
      "dsigma_dt_vm_brodsky_Q2W",
      [](double* xx, double* par) {
        const double t = -std::fabs(xx[0]);
        const double R_n = par[9];
        const double R_c = par[8];
        const double dipole_n = par[7];
        const double Ebeam = par[6];
        const double W = par[5];
        const double s = W * W;
        const double Q2 = par[4];
        const double Mv = par[3];
        const double b = par[2];
        const double c3g = par[1];
        const double c2g = par[0];
        // bail on for unphysical t
        auto tl = pcsim::physics::t_range(s, Q2, M_P, Mv, M_P);
        if (tl.excludes(t)) {
          return 0.;
        }
        return pcsim::physics::dsigma_dt_vm_brodsky(s, t, M_P, Mv, b, c2g,
                                                    c3g) *
               pcsim::physics::dipole_ff_vm(Q2, Mv, dipole_n) *
               (1. + pcsim::physics::R_vm_martynov(Q2, Mv, R_c, R_n) *
                         epsilon(Q2, W, Ebeam));
      },
      (tmin < tmax) ? tmin : tmax, (tmin < tmax) ? tmax : tmin, 10);
}

// =============================================================================
// r00_04 matrix element
// =============================================================================
// x[0]: Q2
// par[0]: W
// par[1]: Mv
// par[2]: Ebeam
// par[3]: R_c
// par[4]: R_n
TF1* r00_04(const double Q2min, const double Q2max) {
  return new TF1("r00_04",
                 [](double* xx, double* par) {
                   const double Q2 = xx[0];
                   const double W = par[0];
                   const double Mv = par[1];
                   const double Ebeam = par[2];
                   const double R_c = par[3];
                   const double R_n = par[4];
                   const double eR =
                       epsilon(Q2, W, Ebeam) *
                       pcsim::physics::R_vm_martynov(Q2, Mv, R_c, R_n);
                   return eR / (1. + eR);
                 },
                 Q2min, Q2max, 5);
}
// =============================================================================
// SCHC theta angle (cos(theta))
// =============================================================================
// x[0]: cos(theta)
// par[0]: Q2
// par[1]: W
// par[2]: Mv
// par[3]: Ebeam
// par[4]: R_c
// par[5]: R_n
TF1* ctheta_schc(const double cthmin = -1, const double cthmax = 1) {
  return new TF1(
      "ctheta_schc",
      [](double* xx, double* par) {
        const double ctheta = xx[0];
        const double Q2 = par[0];
        const double W = par[1];
        const double Mv = par[2];
        const double Ebeam = par[3];
        const double R_c = par[4];
        const double R_n = par[5];
        const double eR = epsilon(Q2, W, Ebeam) *
                          pcsim::physics::R_vm_martynov(Q2, Mv, R_c, R_n);
        const double r04 = eR / (1. + eR);
        return 3. / 8. * ((1. + r04) + (1. - 3. * r04) * ctheta * ctheta);
      },
      cthmin, cthmax, 6);
}
// =============================================================================
// t-integrated VM cross section differential in ctheta as a function of (Q2, W)
// =============================================================================
// x[0]: cos(theta)
// par[0]: 2-gluon Amplitude
// par[1]: 3-gluon Amplitude
// par[2]: b
// par[3]: VM mass
// par[4]: Q2
// par[5]: W (unused)
// par[6]: Ebeam
// par[7]: dipole_n
// par[8]: R_c
// par[9]: R_n
TF1* dsigma_dctheta_vm_brodsky_Q2W(const double cthmin = -1,
                                   const double cthmax = 1) {
  return new TF1(
      "dsigma_dctheta_vm_brodsky_Q2W",
      [=](double* xx, double* par) {
        const double ctheta = xx[0];
        const double R_n = par[9];
        const double R_c = par[8];
        const double dipole_n = par[7];
        const double Ebeam = par[6];
        const double W = par[5];
        const double s = W * W;
        const double Q2 = par[4];
        const double Mv = par[3];
        const double b = par[2];
        const double c3g = par[1];
        const double c2g = par[0];
        // get the t-range, transform to exp(bt)
        auto tl = pcsim::physics::t_range(s, Q2, M_P, Mv, M_P);
        tl.min = std::exp(b * tl.min);
        tl.max = std::exp(b * tl.max);
        // We integrate the dsigma/dexp(bt) formula, which is
        // constant as a function of t. Hence, we only need to
        // evaluate the function once, and multiply with the
        // exp(bt)-range
        const double delta_t = tl.min - tl.max;
        double xs[1] = {tl.min};
        double integral =
            pcsim::physics::dsigma_dexp_bt_vm_brodsky(s, M_P, Mv, b, c2g, c3g) *
            delta_t;
        // add rest of Q2 dependence
        integral *= pcsim::physics::dipole_ff_vm(Q2, Mv, dipole_n) *
                    (1. + pcsim::physics::R_vm_martynov(Q2, Mv, R_c, R_n) *
                              epsilon(Q2, W, Ebeam));
        const double eR = epsilon(Q2, W, Ebeam) *
                          pcsim::physics::R_vm_martynov(Q2, Mv, R_c, R_n);
        const double r04 = eR / (1. + eR);
        const double angular =
            3. / 8. * ((1. + r04) + (1. - 3. * r04) * ctheta * ctheta);
        return integral * angular;
      },
      cthmin, cthmax, 10);
}

} // namespace vm
