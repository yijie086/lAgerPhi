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

#include "fermi.hh"

#include <TMath.h>

extern "C" {
void fermi87_(double& P, int& A, double& res);
}

namespace lager {
namespace physics {

// return 1D fermi distribution as a function of momentum
// in units of GeV^-1
double fermi87(double P, int A) {
  double res;
  double P_MeV = P * 1000;
  fermi87_(P_MeV, A, res);
  res *= 1e9; // MeV^-3 --> GeV^-3
  // Jacobian to go to spherical coordinates, and integrate out theta and phi
  res *= 4 * TMath::Pi() * P * P;
  return res;
}
}
}
