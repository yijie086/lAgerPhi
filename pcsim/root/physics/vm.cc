#include "vm.hh"

// ROOT module with TF1 versions of various physics functions
// see header for full documentation

#include "physics.hh"

#include <algorithm>
#include <cmath>

#include <pcsim/physics/kinematics.hh>
#include <pcsim/physics/photon.hh>
#include <pcsim/physics/vm.hh>

namespace pcsim {
namespace root {
namespace physics {
namespace vm {

TF1* epsilon_Q2(const range_type& Q2lim) {
  return new TF1("epsilon_Q2_WE",
                 [](double* xx, double* par) {
                   const double Q2 = xx[0];
                   const double W = par[0];
                   const double Ebeam = par[1];
                   return epsilon(xx[0], par[0], par[1]);
                 },
                 Q2lim.min, Q2lim.max, 2);
}

TF1* dsigma_dt_vm_brodsky(const range_type& tlim) {
  return new TF1("dsigma_dt_vm_brodsky",
                 [](double* xx, double* par) {
                   const double Mt = M_P;
                   const double t = -std::fabs(xx[0]);
                   const double W = par[4];
                   const double s = W * W;
                   const double Mv = par[3];
                   const double b = par[2];
                   const double c3g = par[1];
                   const double c2g = par[0];
                   return pcsim::physics::dsigma_dt_vm_brodsky(s, t, Mt, Mv, b,
                                                               c2g, c3g);
                 },
                 std::min(tlim.min, tlim.max), std::max(tlim.min, tlim.max), 5);
}

TF1* dsigma_dexp_bt_vm_brodsky(const range_type& tlim) {
  return new TF1("dsigma_dexp_bt_vm_brodsky",
                 [](double* xx, double* par) {
                   const double Mt = M_P;
                   const double t = -std::fabs(xx[0]);
                   const double W = par[4];
                   const double s = W * W;
                   const double Mv = par[3];
                   const double b = par[2];
                   const double c3g = par[1];
                   const double c2g = par[0];
                   return pcsim::physics::dsigma_dexp_bt_vm_brodsky(
                       s, Mt, Mv, b, c2g, c3g);
                 },
                 std::min(tlim.min, tlim.max), std::max(tlim.min, tlim.max), 5);
}
TF1* sigma_vm_brodsky_W(const range_type& Wlim) {
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
                   const double delta_t = tl.max - tl.min;
                   //double xs[1] = {tl.max};
                   const double integral =
                       pcsim::physics::dsigma_dexp_bt_vm_brodsky(s, M_P, Mv, b,
                                                                 c2g, c3g) *
                       delta_t;
                   return integral;
                 },
                 Wlim.min, Wlim.max, 4.);
}
TF1* R_vm_martynov(const range_type& Q2lim) {
  return new TF1("R_vm_martynov",
                 [](double* xx, double* par) {
                   const double Q2 = xx[0];
                   const double Mv = par[0];
                   const double c = par[1];
                   const double n = par[2];
                   return pcsim::physics::R_vm_martynov(Q2, Mv, c, n);
                 },
                 Q2lim.min, Q2lim.max, 3);
}
TF1* dipole_ff_vm(const range_type& Q2lim) {
  return new TF1("dipole_ff_vm",
                 [](double* xx, double* par) {
                   const double Q2 = xx[0];
                   const double Mv = par[0];
                   const double n = par[1];
                   return pcsim::physics::dipole_ff_vm(Q2, Mv, n);
                 },
                 Q2lim.min, Q2lim.max, 2);
}

TF1* sigma_vm_brodsky_W_Q2(const range_type& Wlim) {
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
      Wlim.min, Wlim.max, 10);
}

TF1* dsigma_dt_vm_brodsky_Q2W(const range_type& tlim) {
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
      std::min(tlim.min, tlim.max), std::max(tlim.min, tlim.max), 10);
}

TF1* r00_04(const range_type& Q2lim) {
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
                 Q2lim.min, Q2lim.max, 5);
}
TF1* ctheta_schc(const range_type& cthlim) {
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
      cthlim.min, cthlim.max, 6);
}
TF1* dsigma_dctheta_vm_brodsky_Q2W(const range_type& cthlim) {
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
      cthlim.min, cthlim.max, 10);
}

} // namespace vm
} // namespace physics
} // namespace root
} // namespace pcsim
