// lAger: General Purpose l/A-event Generator
// Copyright (C) 2016-2021 Sylvester Joosten <sjoosten@anl.gov>
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

#include <TMath.h>
#include <cmath>
#include <complex>
#include <gsl/gsl_integration.h>
#include <iostream>

namespace oleksii_jpsi_impl {

//namespace {

constexpr const double kappa = 1.7928473508; // CODATA proton anomalous moment
constexpr const double MD = 1.86484;         // D-meson mass
constexpr const double Mp = 0.938272046;     // proton mass
constexpr const double me = 0.0005109989461; // electron mass
constexpr const double Mj = 3.096916;        // J/psi mass
constexpr const double fj = 0.278;           // J/psi decay constant
constexpr const double Gj = 92.9e-6;         // J/psi total width
constexpr const double h = 4.135667662e-15;  // eV * s (CODATA)
constexpr const double c = 299792458;        // m/s (exact) (CODATA)
const double hc = h * c / (2 * TMath::Pi()) * 1e6; // fm * GeV (CODATA)

constexpr const double alpha = 0.0072973525664; // (CODATA)
const double e = sqrt(4 * TMath::Pi() * alpha); // electric charge
constexpr const double L2 = 0.71; // GeV^2 (SJJ: no idea what this is)

double s2nu(const double s) { return .5 * (s - Mp * Mp - Mj * Mj); }

class Disc {
public:
  constexpr Disc(const double nu0, const double a, const double b,
                 const double C)
      : nu0_{nu0}, a_{a}, b_{b}, C_{C} {}

  double ImT(const double nu) const {
    if (nu < nu0_) {
      return 0;
    }
    return C_ * pow(1. - nu0_ / nu, b_) * pow(nu / nu0_, a_);
  }
  double nu0() const { return nu0_; }

private:
  const double nu0_;
  const double a_;
  const double b_;
  const double C_;
};

struct integrand_param {
  double nu;
  const Disc* el;
  const Disc* inel;
};
double fIntegrand(double nup, void* param) {
  if (!param) {
    return 0.;
  }
  integrand_param* p = static_cast<integrand_param*>(param);
  const double nu = p->nu;
  const Disc* el = p->el;
  const Disc* inel = p->inel;
  const double fnup = el->ImT(nup) + inel->ImT(nup);
  const double fnu = el->ImT(nu) + inel->ImT(nu);
  return (fnup / nup - fnu / nu) / (nup * nup - nu * nu);
}

class jpsi_T_calc {
public:
  jpsi_T_calc()
      : el_{Mp * Mj, 1.38, 1.26, .1}
      , inel_{s2nu(pow(Mp + 2 * MD, 2)), 1.2, 4.2, 18.2}
      , w_{gsl_integration_workspace_alloc(1000)} {
    F_.function = &fIntegrand;
    param_.el = &el_;
    param_.inel = &inel_;
  }

  ~jpsi_T_calc() { gsl_integration_workspace_free(w_); }

  double Im(const double nu) const { return el_.ImT(nu) + inel_.ImT(nu); }
  double Re(const double nu, const double T_0) {
    return T_0 + 2. / TMath::Pi() * nu * nu * disp_int(nu);
  }
  std::complex<double> operator()(const double nu, const double T_0) {
    return {Re(nu, T_0), Im(nu)};
  }

private:
  double disp_int(const double nu) {
    param_.nu = nu;
    double result, error;
    F_.params = &param_;
    gsl_integration_qagiu(&F_, el_.nu0(), 0, 1e-7, 1000, w_, &result, &error);
    return result + F_.function(nu, nullptr) / nu *
                        std::log(fabs((el_.nu0() + nu) / (el_.nu0() - nu))) /
                        (2 * nu);
  }
  const Disc el_;
  const Disc inel_;
  double nu_;
  gsl_integration_workspace* w_;
  integrand_param param_;
  gsl_function F_;
};

jpsi_T_calc jpsi_T;

// Electric form factor
double fGE(const double t) {
  return 1 / std::pow((1 - t / L2), 2);
} // namespace
// Magnetic form factor
double fGM(const double t) { return (1 + kappa) * fGE(t); }
// Pauli form factor
double fF2(const double t) { return kappa * fGE(t) / (1 - t / (4 * Mp * Mp)); }
// Dirac form factor
double fF1(const double t) { return fGM(t) - fF2(t); }


// J/psi resonance T1-amplitude part, as a function of
// q2: invariant mass of dilepton pair squared
// Egamma: photon energy
std::complex<double> fT1(const double q2, const double Egamma, const double t,
                         const double T_0) {
  const std::complex<double> denom = {q2 - Mj * Mj, Mj * Gj};
  const auto fact = fj * fj / (2 * Mp) * 1. / denom;
  const double nu = Mp * Egamma - .5 * Mj * Mj;
  return fact * jpsi_T(nu, T_0) * exp(1.13 * t * 0.5);
}

double Re_jpsi_T(const double q2, const double Egamma, const double T_0) {
  const double nu = Mp * Egamma - .5 * Mj * Mj;
  return std::real(jpsi_T(nu, T_0));
}
double Im_jpsi_T(const double q2, const double Egamma, const double T_0) {
  const double nu = Mp * Egamma - .5 * Mj * Mj;
  return std::imag(jpsi_T(nu, T_0));
}
double Re_jpsi_fT(const double q2, const double Egamma, const double t,
                  const double T_0) {
  return std::real(fT1(q2, Egamma, t, T_0));
}
double Im_jpsi_fT(const double q2, const double Egamma, const double t,
                  const double T_0) {
  return std::imag(fT1(q2, Egamma, t, T_0));
}

// calculate the matrix element squared
// t: mandelstam t
// q2: invariant mass of dilepton pair squared
// Egamma: photon energy
// thetaCM: electron theta in dilepton CM frame
// phiCM: electron phi in dilepton CM frame
double TT(const double t, const double q2, const double Egamma,
          const double thetaCM, const double phiCM, const double T_0) {
  const auto T1 = fT1(q2, Egamma, t, T_0);
  const double ReT1 = std::real(T1);
  const double ImT1 = std::imag(T1);
  const double cth = cos(thetaCM);
  const double sth = sin(thetaCM);
  const double cphi = cos(phiCM);
  const double GEt = fGE(t);
  const double GMt = fGM(t);
  const double Mll = sqrt(q2);
  const double M = Mp;
  const double M2 = M * M;
  const double m = me;
  const double m2 = m * m;
  const double ME = Mp * Egamma;
  const double s = M2 + 2 * ME;
  const double qCM = sqrt((std::pow(2 * ME - q2, 2) - 4 * M2 * q2) / (4 * s));
  const double sqrt_val = sqrt((q2 - 4 * m2) / (M2 + 2 * ME));
  const double a1 = 2 * s * qCM / Mll;
  const double b1 = (q2 * (2 * ME - q2) + t * (2 * ME + q2)) / (2 * Mll * qCM);
  const double b2 = sqrt(s / q2 * std::pow(q2 - t, 2) - b1 * b1);
  const double a = sqrt_val * a1 * cth;
  const double b = sqrt_val * (b1 * cth - b2 * sth * cphi);
  const double L = (Mll * Mll - t - b) * (Mll * Mll - t + b) / 4;
  const double E = Egamma;

  return (std::pow(e, 6) * (std::pow(ImT1, 2) + std::pow(ReT1, 2)) *
          (4 * std::pow(M, 2) - t) *
          (64 * std::pow(M, 4) * std::pow(Mll, 6) +
           4 * std::pow(m, 2) * std::pow(Mll, 8) + std::pow(Mll, 10) -
           128 * std::pow(M, 4) * std::pow(Mll, 4) * t -
           16 * std::pow(m, 2) * std::pow(Mll, 6) * t -
           16 * std::pow(M, 2) * std::pow(Mll, 6) * t -
           4 * std::pow(Mll, 8) * t +
           64 * std::pow(M, 4) * std::pow(Mll, 2) * std::pow(t, 2) +
           24 * std::pow(m, 2) * std::pow(Mll, 4) * std::pow(t, 2) +
           32 * std::pow(M, 2) * std::pow(Mll, 4) * std::pow(t, 2) +
           6 * std::pow(Mll, 6) * std::pow(t, 2) -
           16 * std::pow(m, 2) * std::pow(Mll, 2) * std::pow(t, 3) -
           16 * std::pow(M, 2) * std::pow(Mll, 2) * std::pow(t, 3) -
           4 * std::pow(Mll, 4) * std::pow(t, 3) +
           4 * std::pow(m, 2) * std::pow(t, 4) +
           std::pow(Mll, 2) * std::pow(t, 4) -
           64 * std::pow(m, 2) * M * std::pow(Mll, 6) * E -
           16 * M * std::pow(Mll, 8) * E +
           192 * std::pow(m, 2) * M * std::pow(Mll, 4) * t * E -
           128 * std::pow(M, 3) * std::pow(Mll, 4) * t * E +
           48 * M * std::pow(Mll, 6) * t * E -
           192 * std::pow(m, 2) * M * std::pow(Mll, 2) * std::pow(t, 2) * E +
           128 * std::pow(M, 3) * std::pow(Mll, 2) * std::pow(t, 2) * E -
           16 * M * std::pow(Mll, 4) * std::pow(t, 2) * E +
           64 * std::pow(m, 2) * M * std::pow(t, 3) * E -
           16 * M * std::pow(Mll, 2) * std::pow(t, 3) * E +
           384 * std::pow(m, 2) * std::pow(M, 2) * std::pow(Mll, 4) *
               std::pow(E, 2) +
           96 * std::pow(M, 2) * std::pow(Mll, 6) * std::pow(E, 2) -
           768 * std::pow(m, 2) * std::pow(M, 2) * std::pow(Mll, 2) * t *
               std::pow(E, 2) +
           256 * std::pow(M, 4) * std::pow(Mll, 2) * t * std::pow(E, 2) -
           192 * std::pow(M, 2) * std::pow(Mll, 4) * t * std::pow(E, 2) +
           384 * std::pow(m, 2) * std::pow(M, 2) * std::pow(t, 2) *
               std::pow(E, 2) +
           32 * std::pow(M, 2) * std::pow(Mll, 2) * std::pow(t, 2) *
               std::pow(E, 2) -
           1024 * std::pow(m, 2) * std::pow(M, 3) * std::pow(Mll, 2) *
               std::pow(E, 3) -
           256 * std::pow(M, 3) * std::pow(Mll, 4) * std::pow(E, 3) +
           1024 * std::pow(m, 2) * std::pow(M, 3) * t * std::pow(E, 3) +
           256 * std::pow(M, 3) * std::pow(Mll, 2) * t * std::pow(E, 3) +
           1024 * std::pow(m, 2) * std::pow(M, 4) * std::pow(E, 4) +
           256 * std::pow(M, 4) * std::pow(Mll, 2) * std::pow(E, 4) +
           std::pow(b, 2) * (std::pow(Mll, 6) +
                             16 * std::pow(M, 2) * (4 * std::pow(M, 2) - t) *
                                 std::pow(E, 2) -
                             2 * std::pow(Mll, 4) * (t + 4 * M * E) +
                             std::pow(Mll, 2) * std::pow(t + 4 * M * E, 2)) +
           4 * a * b *
               (std::pow(Mll, 6) + 4 * M * t * (-4 * std::pow(M, 2) + t) * E -
                2 * std::pow(Mll, 4) * (t + 4 * M * E) +
                std::pow(Mll, 2) * (std::pow(t, 2) + 4 * M * t * E +
                                    16 * std::pow(M, 2) * E * (M + E))) +
           4 * std::pow(a, 2) *
               (std::pow(std::pow(Mll, 2) - t, 3) +
                8 * M * std::pow(Mll, 2) * (-std::pow(Mll, 2) + t) * E +
                4 * std::pow(M, 2) *
                    (std::pow(Mll, 4) + std::pow(t, 2) -
                     2 * std::pow(Mll, 2) * (t - 2 * std::pow(E, 2)))))) /
         (std::pow(Mll, 4) * std::pow(-std::pow(Mll, 2) + t + 4 * M * E, 4));
}
//} // namespace

// calculate the differential cross section
// in nb/GeV^4/sr
// t: mandelstam t
// q2: invariant mass of dilepton pair squared
// Egamma: photon energy
// thetaCM: electron theta in dilepton CM frame
// phiCM: electron phi in dilepton CM frame

double calc_xsec(const double t, const double q2, const double Egamma,
                 const double thetaCM, const double phiCM, const double T_0) {

  // TODO NEED TO DOUBLE CHECK THE UNITS AND UNIT CONVERSION TODO
  const double me2 = me * me;
  const double xsec_mb_per_GeV4_sr =
      std::sqrt(1 - 4 * me2 / q2) /
      (1024. * std::pow(TMath::Pi(), 4) * std::pow(2 * Mp * Egamma, 2)) *
      TT(t, q2, Egamma, thetaCM, phiCM, T_0) * hc * hc * 1e4;
  return xsec_mb_per_GeV4_sr * 1e3;
}

} // namespace oleksii_jpsi_impl

