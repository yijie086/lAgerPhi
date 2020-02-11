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

#include <TMath.h>
#include <cmath>
#include <iostream>
#include <complex>
#include <gsl/gsl_integration.h>

namespace oleksii_bh_impl {

namespace {

constexpr const double kappa = 1.7928473508; // CODATA proton anomalous moment
constexpr const double MD = 1.86484;         // D-meson mass
constexpr const double Mp = 0.938272046;     // proton mass
constexpr const double me = 0.0005109989461; // electron mass
constexpr const double Mj = 3.096916;        // J/psi mass
constexpr const double fj = 0.278;           // J/psi decay constant
constexpr const double Gj = 92.9e-6;         // J/psi total width
constexpr const double T_0 = 0.;             // T(0) subtraction constant
constexpr const double h = 4.135667662e-15;  // eV * s (CODATA)
constexpr const double c = 299792458	;// m/s (exact) (CODATA)
const double hc = h * c / (2 * TMath::Pi()) * 1e6; // fm * GeV (CODATA)

constexpr const double alpha = 0.0072973525664; // (CODATA)
const double e = sqrt(4 * TMath::Pi() * alpha); // electric charge
constexpr const double L2 = 0.71; // GeV^2 (SJJ: no idea what this is)

double s2nu(const double s) { return .5 * (s - Mp * Mp - Mj * Mj); }

// Electric form factor
double fGE(const double t) {
  return 1/std::pow((1 - t/L2), 2);
}
// Magnetic form factor
double fGM(const double t) {
  return (1+kappa)*fGE(t);
}
// Pauli form factor
double fF2(const double t) { return kappa * fGE(t) / (1 - t / (4 * Mp * Mp)); }
// Dirac form factor
double fF1(const double t) { return fGM(t) - fF2(t); }

// calculate the matrix element squared
// t: mandelstam t
// q2: invariant mass of dilepton pair squared
// Egamma: photon energy
// thetaCM: electron theta in dilepton CM frame
// phiCM: electron phi in dilepton CM frame
double BH(const double t, const double q2, const double Egamma,
                 const double thetaCM, const double phiCM) {
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

  return (std::pow(e,6)*(-256*std::pow(E,2)*std::pow(GEt,2)*L*std::pow(m,2)*std::pow(M,4) + 
       128*E*std::pow(GEt,2)*L*std::pow(m,2)*std::pow(M,3)*std::pow(Mll,2) - 
       128*std::pow(E,2)*std::pow(GEt,2)*L*std::pow(M,4)*std::pow(Mll,2) - 
       16*std::pow(GEt,2)*L*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,4) + 
       16*std::pow(GMt,2)*L*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,4) + 
       64*E*std::pow(GEt,2)*L*std::pow(M,3)*std::pow(Mll,4) + 
       64*std::pow(E,2)*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,4)*std::pow(Mll,4) - 
       8*std::pow(GEt,2)*L*std::pow(M,2)*std::pow(Mll,6) + 
       8*std::pow(GMt,2)*L*std::pow(M,2)*std::pow(Mll,6) - 
       32*E*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,3)*std::pow(Mll,6) + 
       32*std::pow(E,2)*std::pow(GEt,2)*std::pow(M,4)*std::pow(Mll,6) + 
       4*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,8) - 
       4*std::pow(GMt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,8) - 
       16*E*std::pow(GEt,2)*std::pow(M,3)*std::pow(Mll,8) - 
       8*std::pow(GEt,2)*std::pow(M,4)*std::pow(Mll,8) + 
       2*std::pow(GEt,2)*std::pow(M,2)*std::pow(Mll,10) - 
       2*std::pow(GMt,2)*std::pow(M,2)*std::pow(Mll,10) + 
       64*std::pow(E,2)*std::pow(GMt,2)*L*std::pow(m,2)*std::pow(M,2)*t + 
       256*std::pow(GMt,2)*L*std::pow(m,4)*std::pow(M,2)*t - 
       128*E*std::pow(GEt,2)*L*std::pow(m,2)*std::pow(M,3)*t + 
       256*std::pow(GEt,2)*L*std::pow(m,2)*std::pow(M,4)*t - 
       32*E*std::pow(GMt,2)*L*std::pow(m,2)*M*std::pow(Mll,2)*t + 
       32*std::pow(E,2)*std::pow(GMt,2)*L*std::pow(M,2)*std::pow(Mll,2)*t + 
       32*std::pow(GEt,2)*L*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,2)*t - 
       96*std::pow(GMt,2)*L*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,2)*t - 
       64*E*std::pow(GEt,2)*L*std::pow(M,3)*std::pow(Mll,2)*t - 
       128*std::pow(GEt,2)*L*std::pow(M,4)*std::pow(Mll,2)*t - 
       128*std::pow(E,2)*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,4)*std::pow(Mll,2)*t - 
       16*E*std::pow(GMt,2)*L*M*std::pow(Mll,4)*t + 
       16*std::pow(GEt,2)*L*std::pow(M,2)*std::pow(Mll,4)*t - 
       16*std::pow(GMt,2)*L*std::pow(M,2)*std::pow(Mll,4)*t - 
       16*std::pow(E,2)*std::pow(GMt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,4)*t + 
       64*std::pow(GMt,2)*std::pow(m,4)*std::pow(M,2)*std::pow(Mll,4)*t + 
       96*E*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,3)*std::pow(Mll,4)*t - 
       96*std::pow(E,2)*std::pow(GEt,2)*std::pow(M,4)*std::pow(Mll,4)*t + 
       64*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,4)*std::pow(Mll,4)*t + 
       8*E*std::pow(GMt,2)*std::pow(m,2)*M*std::pow(Mll,6)*t - 
       8*std::pow(E,2)*std::pow(GMt,2)*std::pow(M,2)*std::pow(Mll,6)*t - 
       16*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,6)*t + 
       64*E*std::pow(GEt,2)*std::pow(M,3)*std::pow(Mll,6)*t + 
       32*std::pow(GEt,2)*std::pow(M,4)*std::pow(Mll,6)*t + 4*E*std::pow(GMt,2)*M*std::pow(Mll,8)*t - 
       8*std::pow(GEt,2)*std::pow(M,2)*std::pow(Mll,8)*t + 
       6*std::pow(GMt,2)*std::pow(M,2)*std::pow(Mll,8)*t - 
       64*std::pow(GMt,2)*L*std::pow(m,4)*std::pow(t,2) + 
       32*E*std::pow(GMt,2)*L*std::pow(m,2)*M*std::pow(t,2) - 
       80*std::pow(GEt,2)*L*std::pow(m,2)*std::pow(M,2)*std::pow(t,2) + 
       16*std::pow(GMt,2)*L*std::pow(m,2)*std::pow(M,2)*std::pow(t,2) + 
       64*std::pow(E,2)*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,4)*std::pow(t,2) + 
       16*std::pow(GMt,2)*L*std::pow(m,2)*std::pow(Mll,2)*std::pow(t,2) + 
       16*E*std::pow(GMt,2)*L*M*std::pow(Mll,2)*std::pow(t,2) + 
       24*std::pow(GEt,2)*L*std::pow(M,2)*std::pow(Mll,2)*std::pow(t,2) - 
       24*std::pow(GMt,2)*L*std::pow(M,2)*std::pow(Mll,2)*std::pow(t,2) + 
       32*std::pow(E,2)*std::pow(GMt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,2)*std::pow(t,2) - 
       128*std::pow(GMt,2)*std::pow(m,4)*std::pow(M,2)*std::pow(Mll,2)*std::pow(t,2) - 
       96*E*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,3)*std::pow(Mll,2)*std::pow(t,2) + 
       96*std::pow(E,2)*std::pow(GEt,2)*std::pow(M,4)*std::pow(Mll,2)*std::pow(t,2) - 
       128*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,4)*std::pow(Mll,2)*std::pow(t,2) - 
       16*std::pow(GMt,2)*std::pow(m,4)*std::pow(Mll,4)*std::pow(t,2) - 
       24*E*std::pow(GMt,2)*std::pow(m,2)*M*std::pow(Mll,4)*std::pow(t,2) + 
       24*std::pow(E,2)*std::pow(GMt,2)*std::pow(M,2)*std::pow(Mll,4)*std::pow(t,2) + 
       8*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,4)*std::pow(t,2) + 
       40*std::pow(GMt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,4)*std::pow(t,2) - 
       96*E*std::pow(GEt,2)*std::pow(M,3)*std::pow(Mll,4)*std::pow(t,2) - 
       48*std::pow(GEt,2)*std::pow(M,4)*std::pow(Mll,4)*std::pow(t,2) + 
       4*std::pow(GMt,2)*std::pow(m,2)*std::pow(Mll,6)*std::pow(t,2) - 
       16*E*std::pow(GMt,2)*M*std::pow(Mll,6)*std::pow(t,2) + 
       12*std::pow(GEt,2)*std::pow(M,2)*std::pow(Mll,6)*std::pow(t,2) - 
       4*std::pow(GMt,2)*std::pow(M,2)*std::pow(Mll,6)*std::pow(t,2) + 
       std::pow(GMt,2)*std::pow(Mll,8)*std::pow(t,2) - 
       16*std::pow(E,2)*std::pow(GMt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(t,3) + 
       64*std::pow(GMt,2)*std::pow(m,4)*std::pow(M,2)*std::pow(t,3) + 
       32*E*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,3)*std::pow(t,3) - 
       32*std::pow(E,2)*std::pow(GEt,2)*std::pow(M,4)*std::pow(t,3) + 
       64*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,4)*std::pow(t,3) + 
       8*std::pow(GMt,2)*L*std::pow(Mll,2)*std::pow(t,3) + 
       32*std::pow(GMt,2)*std::pow(m,4)*std::pow(Mll,2)*std::pow(t,3) + 
       24*E*std::pow(GMt,2)*std::pow(m,2)*M*std::pow(Mll,2)*std::pow(t,3) - 
       24*std::pow(E,2)*std::pow(GMt,2)*std::pow(M,2)*std::pow(Mll,2)*std::pow(t,3) + 
       16*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,2)*std::pow(t,3) - 
       64*std::pow(GMt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(Mll,2)*std::pow(t,3) + 
       64*E*std::pow(GEt,2)*std::pow(M,3)*std::pow(Mll,2)*std::pow(t,3) + 
       32*std::pow(GEt,2)*std::pow(M,4)*std::pow(Mll,2)*std::pow(t,3) - 
       16*std::pow(GMt,2)*std::pow(m,2)*std::pow(Mll,4)*std::pow(t,3) + 
       24*E*std::pow(GMt,2)*M*std::pow(Mll,4)*std::pow(t,3) - 
       8*std::pow(GEt,2)*std::pow(M,2)*std::pow(Mll,4)*std::pow(t,3) - 
       4*std::pow(GMt,2)*std::pow(M,2)*std::pow(Mll,4)*std::pow(t,3) - 
       4*std::pow(GMt,2)*std::pow(Mll,6)*std::pow(t,3) - 16*std::pow(GMt,2)*std::pow(m,4)*std::pow(t,4) - 
       8*E*std::pow(GMt,2)*std::pow(m,2)*M*std::pow(t,4) + 
       8*std::pow(E,2)*std::pow(GMt,2)*std::pow(M,2)*std::pow(t,4) - 
       12*std::pow(GEt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(t,4) + 
       28*std::pow(GMt,2)*std::pow(m,2)*std::pow(M,2)*std::pow(t,4) - 
       16*E*std::pow(GEt,2)*std::pow(M,3)*std::pow(t,4) - 8*std::pow(GEt,2)*std::pow(M,4)*std::pow(t,4) + 
       20*std::pow(GMt,2)*std::pow(m,2)*std::pow(Mll,2)*std::pow(t,4) - 
       16*E*std::pow(GMt,2)*M*std::pow(Mll,2)*std::pow(t,4) + 
       2*std::pow(GEt,2)*std::pow(M,2)*std::pow(Mll,2)*std::pow(t,4) + 
       6*std::pow(GMt,2)*std::pow(M,2)*std::pow(Mll,2)*std::pow(t,4) + 
       6*std::pow(GMt,2)*std::pow(Mll,4)*std::pow(t,4) - 8*std::pow(GMt,2)*std::pow(m,2)*std::pow(t,5) + 
       4*E*std::pow(GMt,2)*M*std::pow(t,5) - 2*std::pow(GMt,2)*std::pow(M,2)*std::pow(t,5) - 
       4*std::pow(GMt,2)*std::pow(Mll,2)*std::pow(t,5) + std::pow(GMt,2)*std::pow(t,6) + 
       4*std::pow(a,2)*(std::pow(b,2)*std::pow(m,2) + L*(4*std::pow(m,2) - 2*t) + 
          std::pow(m,2)*std::pow(std::pow(Mll,2) - t,2))*
        (4*std::pow(GEt,2)*std::pow(M,2) - std::pow(GMt,2)*t) + 
       a*b*(4*std::pow(GEt,2)*std::pow(M,2) - std::pow(GMt,2)*t)*
        (4*L*(4*std::pow(m,2) + 4*E*M - std::pow(Mll,2) - t) + 
          (-4*E*M*(8*std::pow(m,2) + std::pow(Mll,2) - t) + 
             (std::pow(Mll,2) - t)*(12*std::pow(m,2) + std::pow(Mll,2) - t))*(std::pow(Mll,2) - t) + 
          std::pow(b,2)*(4*std::pow(m,2) + 4*E*M - std::pow(Mll,2) + t)) + 
       std::pow(b,4)*(2*std::pow(GEt,2)*std::pow(M,2)*
           (2*std::pow(m,2) + 4*E*M + 4*std::pow(M,2) - std::pow(Mll,2)) - 
          std::pow(GMt,2)*(4*std::pow(m,2)*std::pow(M,2) + 2*std::pow(M,2)*(std::pow(Mll,2) - t) + 
             2*E*M*t + t*(-std::pow(Mll,2) + t))) + 
       std::pow(b,2)*(8*std::pow(GEt,2)*std::pow(M,2)*
           (L*(2*std::pow(m,2) + 4*E*M - std::pow(Mll,2)) + 
             E*M*(std::pow(Mll,2) - t)*(-12*std::pow(m,2) + std::pow(Mll,2) - t) + 
             4*std::pow(E,2)*std::pow(M,2)*(2*std::pow(m,2) - std::pow(Mll,2) + t) + 
             std::pow(m,2)*(3*std::pow(Mll,4) + 8*std::pow(M,2)*t - 6*std::pow(Mll,2)*t + std::pow(t,2))
             ) - std::pow(GMt,2)*(16*std::pow(m,4)*t*(-4*std::pow(M,2) + t) + 
             4*L*(4*std::pow(m,2)*std::pow(M,2) + 2*std::pow(M,2)*std::pow(Mll,2) + 2*E*M*t - 
                std::pow(Mll,2)*t) + (std::pow(Mll,2) - t)*
              (2*E*M*(std::pow(Mll,2) - t)*t + std::pow(Mll,2)*(std::pow(Mll,2) - t)*t - 
                4*std::pow(M,2)*(std::pow(Mll,4) + 2*std::pow(E,2)*t - std::pow(Mll,2)*t)) - 
             4*std::pow(m,2)*(6*E*M*(std::pow(Mll,2) - t)*t + 
                t*(-2*std::pow(Mll,4) + 3*std::pow(Mll,2)*t - 2*std::pow(t,2)) + 
                2*std::pow(M,2)*(std::pow(Mll,4) - 2*std::pow(E,2)*t + std::pow(t,2)))))))/
   (2.*std::pow(L,2)*(4*std::pow(M,2) - t)*std::pow(t,2));

}
} // namespace

// calculate the differential cross section
// in nb/GeV^4/sr
// t: mandelstam t
// q2: invariant mass of dilepton pair squared
// Egamma: photon energy
// thetaCM: electron theta in dilepton CM frame
// phiCM: electron phi in dilepton CM frame


double calc_xsec(const double t, const double q2, const double Egamma,
                 const double thetaCM, const double phiCM) {

  // TODO NEED TO DOUBLE CHECK THE UNITS AND UNIT CONVERSION TODO
  const double me2 = me * me;
  const double xsec_mb_per_GeV4_sr =
      std::sqrt(1 - 4 * me2 / q2) /
      (1024. * std::pow(TMath::Pi(), 4) * std::pow(2 * Mp * Egamma, 2)) *
      BH(t, q2, Egamma, thetaCM, phiCM) * hc * hc * 1e4;
  return xsec_mb_per_GeV4_sr * 1e3;
}

} // namespace oleksii_bh_impl
