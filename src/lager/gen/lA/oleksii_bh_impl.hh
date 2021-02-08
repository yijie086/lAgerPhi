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

#ifndef OLEKSII_BH_IMPL_LOADED
#define OLEKSII_BH_IMPL_LOADED

namespace oleksii_bh_impl {

// calculate the cross section
// t: mandelstam t
// q2: invariant mass of dilepton pair squared
// Egamma: photon energy
// thetaCM: electron theta in dilepton CM frame
// phiCM: electron phi in dilepton CM frame
double calc_xsec(const double t, const double q2, const double Egamma,
                 const double thetaCM, const double phiCM);
} // namespace oleksii_bh_impl

#endif
