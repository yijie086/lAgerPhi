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

#ifndef OLEKSII_JPSI_IMPL_LOADED
#define OLEKSII_JPSI_IMPL_LOADED

namespace oleksii_jpsi_impl {

std::complex<double> fT1(const double q2, const double Egamma);

double TT(const double t, const double q2, const double Egamma,
          const double thetaCM, const double phiCM);

double Re_jpsi_T(const double q2, const double Egamma, const double T_0 = 0);
double Im_jpsi_T(const double q2, const double Egamma, const double T_0 = 0);
double Re_jpsi_fT(const double q2, const double Egamma, const double t,
                  const double T_0 = 0);
double Im_jpsi_fT(const double q2, const double Egamma, const double t,
                  const double T_0 = 0);

// calculate the cross section
// t: mandelstam t
// q2: invariant mass of dilepton pair squared
// Egamma: photon energy
// thetaCM: electron theta in dilepton CM frame
// phiCM: electron phi in dilepton CM frame
double calc_xsec(const double t, const double q2, const double Egamma,
                 const double thetaCM, const double phiCM,
                 const double T_0 = 0);
} // namespace oleksii_jpsi_impl

#endif
